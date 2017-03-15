#!/usr/bin/env python2

###########################################################
# Automatic mapping of WGS reads to a library of RAD loci #
###########################################################


import os
import sys
import glob
import subprocess
import gzip
import re
import operator
from optparse import OptionParser
import multiprocessing


def create_shotgun2rad_params(version):
    """
    Creates a new shotgun2rad_params.txt file
    """

    output = """
==== parameter inputs for shotgun2rad version %s  ===========================================
pyRAD                      ## 1.: command (or path) to call pyRAD
bwa                        ## 2.: command (or path) to call bwa
1                          ## 3.: number of threads for bwa (def. 1)
samtools                   ## 4.: command (or path) to call samtools
                           ## 5.: path and prefix of bwa-indexed loci of reference
0.90                       ## 6.: minOverlap: of read mapped against ref. locus (def. 0.9)
60                         ## 7.: minMAPQ: minimum mapping quality (0-60, def. 60)
                           ## 8.: clustSize: reads to retain per cluster (def. = mindepth)
==========================================================================================""" % (version)
    outfile = open("shotgun2rad_params.txt", 'w')
    print >> sys.stderr, "\tnew wg2rad_params.txt file created"
    print >> outfile, "\n".join(output.split("\n")[1:])


def cmd_exists(cmd):
    """
    Checks existence of a program

    :param cmd: name or path/name of the program
    """
    return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0


def mapeditedreads(editfile, bwa, bwa_threads, mapref, samtools):
    """
    Maps quality-filtered reads to loci of reference

    :param editfile: name of reads file edited by pyRAD's step 2
    :param bwa: command or path for bwa
    :param bwa_threads: number of threads for bwa, will divide number of processors in pyRAD's params
    :param mapref: location and prefix of bwa-indexed reference
    :param samtools: command or path for samtools 1.3.1
    :return: mapped reads in SAM format in WORK/maps/ folder
    """

    samplein = editfile.replace(".edit", "")
    mapsample = samplein.replace("/edits", "/maps")
    map_cmd = bwa + " mem -t " + str(
        bwa_threads) + " " + mapref + " " + samplein + ".edit | " + samtools + " view -hS -F 4 - >" + mapsample + ".step"
    subprocess.call(map_cmd, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    sort_cmd = samtools + " sort -o " + mapsample + ".sam -O sam -@ " + str(bwa_threads) + " " + mapsample + ".step"
    subprocess.call(sort_cmd, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    os.remove(mapsample + ".step")
    sys.stderr.write(".")


def parsecigar(cigarstring, seq):
    """
    Modifies a sequence according to its CIGAR string

    :param cigarstring: cigar string, check SAM format specification
    :param seq: raw sequence
    :return: edited sequence according to cigar string
    """

    matches = re.findall(r'(\d+)([A-Z]{1})', cigarstring)
    cigar = [{'type': m[1], 'length': int(m[0])} for m in matches]
    start = 0
    seqout = ''
    for i in range(0, len(cigar)):
        l = cigar[i]['length']
        if cigar[i]['type'] == 'S':
            start += l
        elif cigar[i]['type'] == 'H':
            continue
        elif cigar[i]['type'] == 'M':
            seqout += seq[start:start + l]
            start += l
        elif cigar[i]['type'] == 'I':
            seqout += seq[start:start + l]
            start += l
        elif cigar[i]['type'] == 'D':
            seqout += '-' * l
        else:
            print "SAM file probably contains unmapped reads"

    return seqout


def alignfast(names, seqs, muscle):
    """
    Performs muscle alignments on cluster and returns output as string,
    function taken directly from pyRAD 3.0.5

    :param names: list with names of sequences
    :param seqs: list with sequences
    :param muscle: muscle alias or path
    :return: raw muscle alignment in a single string
    """

    ST = "\n".join('>' + i + '\n' + j for i, j in zip(names, seqs))
    cmd = "/bin/echo '" + ST + "' | " + muscle + " -quiet -in -"
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (fin, fout) = (p.stdin, p.stdout)

    return fout.read()


def parsemuscle(musclestdout):
    """
    Parses muscle output from a string to two lists,
    function taken directly from pyRAD 3.0.5 [ sortalign() ]

    :param musclestdout: raw muscle stdout as a single string
    :return: list of sequence names, list of sequences
    """

    G = musclestdout.split("\n>")
    GG = [i.split("\n")[0].replace(">", "") + "\n" + "".join(i.split("\n")[1:]) for i in G]
    alignment = [i.split("\n") for i in GG]
    names = [">" + i[0] for i in alignment]
    seqs = [i[1] for i in alignment]

    return names, seqs


def trimalignment(untrimmed_alignment):
    """
    Trims leading and trailing "-" from an alignment

    :param untrimmed_alignment: untrimmed alignment
    :return: trimmed alignment
    """

    maxtrimright = 0
    maxtrimleft = 0
    aln_len = len(untrimmed_alignment[0])
    for seq in untrimmed_alignment:
        if len(seq) - len(seq.lstrip("-")) > maxtrimleft:
            maxtrimleft = len(seq) - len(seq.lstrip("-"))
        if len(seq) - len(seq.rstrip("-")) > maxtrimright:
            maxtrimright = len(seq) - len(seq.rstrip("-"))
    trimmed_alignment = [seq[(maxtrimleft//4):aln_len - (maxtrimright//4)] for seq in untrimmed_alignment]

    return trimmed_alignment


def derepsortbysize(names, seqs, clustsize):
    """
    Dereplicates, counts duplicates, and sorts alignment based
    on the abundance of the sequence. Number of sequences returned
    is at least clustsize

    :param names: list of sequence names
    :param seqs: list of sequences
    :return: alignment sorted by abundance as two lists:
             list of sequences names with abundance in ";size=n;"
    """

    derep = {}
    for i, j in zip(names, seqs):
        if j in derep.keys():
            derep[j][0] = i
            derep[j][1] += 1
        else:
            derep[j] = [i, 1]

    names = []
    seqs = []
    sizes = []
    for seq in derep:
        names.append(derep[seq][0] + ";size=" + str(derep[seq][1]) + ";")
        seqs.append(seq)
        sizes.append(derep[seq][1])

    sortedbysize_full = sorted(zip(names, seqs, sizes), reverse=True, key=lambda x: x[2])

    sortedbysize = []
    depth = 0
    for k in range(len(sortedbysize_full)):
        sortedbysize.append(sortedbysize_full[k])
        depth += sortedbysize_full[k][2]
        if depth >= clustsize:
            break

    read_names, read_seqs = [[x[l] for x in sortedbysize] for l in range(2)]

    return read_names, read_seqs


def removeholes(seqs):
    """
    Removes inner deletions shared by all sequences in an alignment

    :param seqs: alignment with empty columns
    :return: alignment without empty columns
    """

    consensus = list('-' * len(seqs[0]))
    for seq in seqs:
        for i in range(len(seq)):
            if seq[i] != '-':
                consensus[i] = 'N'

    positions = []
    for j in range(len(consensus)):
        if consensus[j] == '-':
            positions.append(j)
    positions = sorted(positions, reverse=True)

    newseqs = []
    for seq in seqs:
        newseq = seq
        for pos in positions:
            newseq = newseq[:pos] + newseq[pos + 1:]
        newseqs.append(newseq)

    return newseqs


def makeclustS(mapsam, outfolder, minMAPQ, minoverlap, muscle, clustsize):
    """
    Parses SAM files and extract alignments to make clusters file compatible with pyRAD

    :param mapsam: name of SAM file produced by step 2
    :param outfolder: destination of clusters, WORK/clust.XX
    :param minMAPQ: minimum mapping quality for each read
    :param minoverlap: minimum overlap proportion of read mapped against reference locus
    :param muscle: command or path for muscle
    :param clustsize: minimal sequencing depth to form a cluster
    :return: sample.clustS.gz file in outfolder
    """

    with open(mapsam) as sam:
        clustfilename = outfolder + "/" + mapsam.split("/")[-1].replace(".sam", ".clustS.gz")
        outclust = gzip.open(clustfilename, 'w')
        previouslocus = ''
        clusters = []
        loci_lengths = {}
        read_names = []
        read_seqs = []
        map_quals = []
        insertions = 0
        deletions = 0
        for line in sam:

            """ store reference loci names and their lengths in a dictionary """
            if line[0:3] == "@SQ":
                sn = line.strip("\n").split("\t")[1].replace("SN:", "")
                ln = int(line.strip("\n").split("\t")[2].replace("LN:", ""))
                loci_lengths[sn] = ln

                """ verify mapping quality according to minMAPQ parameter """
            elif line[0] != "@":
                mapq = int(line.split("\t")[4])
                if mapq >= minMAPQ:
                    readname = line.split("\t")[0]
                    locusname = line.split("\t")[2]
                    position = int(line.split("\t")[3])
                    cigar = line.split("\t")[5]
                    seqin = line.split("\t")[9]

                    seqout = parsecigar(cigar, seqin)
                    trailing = loci_lengths[locusname] - (
                    len(seqout) + position - 1)  # number of trailing "-" to pad sequence

                    """ check if the cigar-edited sequence covers reference in at least minOverlap proportion """
                    if len(seqout) >= loci_lengths[locusname] * minoverlap:
                        seqout = '-' * (
                        position - 1) + seqout + '-' * trailing  # pad aligned sequence with "-" on both ends

                        """ check if still extracting reads mapped to the same locus """
                        if locusname == previouslocus or previouslocus == '':
                            if 'I' in cigar:
                                insertions += 1
                            if 'D' in cigar:
                                deletions += 1
                            read_names.append(readname)
                            read_seqs.append(seqout)
                            map_quals.append(mapq)

                            """ or if it is time to append alignment to list of clusters """
                        else:

                            """ when alignment consists of a single sequence """
                            if len(read_seqs) == 1:
                                clusters.append("\n".join(">" + i + ";size=1;\n" + j.replace("-", "") for i, j in
                                                          zip(read_names, read_seqs)))
                                read_names = [readname]  # first read of next locus alignment
                                read_seqs = [seqout]  # first read of next locus alignment
                                map_quals = [mapq]  # first read of next locus alignment
                                insertions = 0
                                deletions = 0
                                if 'I' in cigar:
                                    insertions += 1
                                if 'D' in cigar:
                                    deletions += 1

                                """ when alignment has more than one sequence """
                            elif len(read_seqs) > 1:

                                """ sort mapped and aligned reads by mapping quality """
                                sortedbyquality = sorted(zip(map_quals, read_names, read_seqs), reverse=True,
                                                         key=lambda x: x[0])
                                map_quals, read_names, read_seqs = [[x[i] for x in sortedbyquality] for i in
                                                                    range(3)]

                                """ if sequences contained insertions or alignment
                                    had too many deletions will be realigned with muscle"""
                                if insertions > 0 or deletions > len(read_seqs) * 0.25:

                                    """ if alignment has more than 200 sequences, send only the set of
                                        200 sequences with best mapping quality to muscle """
                                    if len(read_seqs) > 500:
                                        read_names, read_seqs = read_names[0:500], read_seqs[0:500]
                                    muscle_string = alignfast(read_names, read_seqs, muscle)
                                    read_names, read_seqs = parsemuscle(muscle_string)
                                    read_names, read_seqs = derepsortbysize(read_names, read_seqs,
                                                                            clustsize)
                                    read_seqs = removeholes(read_seqs)
                                    clusters.append("\n".join(i + "\n" + j for i, j in zip(read_names, read_seqs)))

                                    """ or if it is a good enough alignment to skip muscle """
                                else:
                                    read_names, read_seqs = derepsortbysize(read_names, read_seqs,
                                                                            clustsize)
                                    read_seqs = removeholes(read_seqs)
                                    clusters.append(
                                        "\n".join(">" + i + "\n" + j for i, j in zip(read_names, read_seqs)))
                                read_names = [readname]  # first read of next locus alignment
                                read_seqs = [seqout]  # first read of next locus alignment
                                map_quals = [mapq]  # first read of next locus alignment
                                insertions = 0
                                deletions = 0
                                if 'I' in cigar:
                                    insertions += 1
                                if 'D' in cigar:
                                    deletions += 1
                        previouslocus = locusname
        outclust.write("\n//\n//\n".join(clusters) + "\n//\n//\n")
        outclust.close()
        sys.stderr.write(".")


class Worker(multiprocessing.Process):
    def __init__(self, work_queue, result_queue, func):

        # base class initialization
        multiprocessing.Process.__init__(self)

        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.func = func

    def run(self):
        while not self.kill_received:
            # get a task
            if self.work_queue.empty():
                break
            else:
                job = self.work_queue.get()

            # the actual processing
            res = self.func(*job)

            # store the result
            self.result_queue.put(res)


def main():
    parser = OptionParser(prog="shotgun2rad", usage="%prog [options]", version="%prog 2.0")
    parser.add_option('-p', action="store", type="string", dest="pyrad_params",
                      help="pyRAD parameters file\n")
    parser.add_option('-w', action="store", type="string", dest="shotgun2rad_params",
                      help="shotgun2rad parameters file\n")
    parser.add_option('-s', action="store", dest="steps",
                      help="""perform step-wise parts of within analysis\n
                      1 = filter/edit raw sequences   (pyrad -s 2) \
                      2 = mapping reads to reference loci          \
                      3 = parsing maps to create clusters          \
                      4 = estimate pi and e           (pyrad -s 4) \
                      5 = consensus calling           (pyrad -s 5) \
                      6 = cluster consensus           (pyrad -s 6) \
                      7 = align & create output files (pyrad -s 7)""")
    parser.add_option('-n', action="store_true", dest="new_shotgun2radparams_file",
                      help="""creates a new empty input shotgun2rad_params.txt file """)

    (options, args) = parser.parse_args()

    if not any([options.pyrad_params, options.shotgun2rad_params, options.new_shotgun2radparams_file]):
        print "\n\tmust include option of -p and -w, or -n (use -h for help)\n"
        sys.exit()

    if options.pyrad_params and options.shotgun2rad_params:
        sys.stderr.write('\n\n' + ' ' * 5 + '---' * 23 + '\n' + \
                         ' ' * 6 + 'shotgun2rad : Whole-Genome Shotgun data to RAD-Seq converter for pyRAD\n' + \
                         ' ' * 5 + '---' * 23 + '\n\n')

        params_pyrad = str(options.pyrad_params)

        """ read shared parameters from pyRAD's params.txt file """
        read_pyrad_params = [line.strip().split('##')[0].strip() for line in open(options.pyrad_params).readlines()]
        if "==** " not in str(read_pyrad_params[0]):
            print "\n\twarning: update params input file format to latest version\n";
            sys.exit()

        WORK     = str(read_pyrad_params[1])
        muscle = str(read_pyrad_params[5])
        parallel = int(read_pyrad_params[7])
        mindepth = int(read_pyrad_params[8])
        wclust   = str(read_pyrad_params[10])

        """ read parameters exclusive for shotgun2rad """
        read_shotgun2rad_params = [line.strip().split('##')[0].strip() for line in open(options.shotgun2rad_params).readlines()]

        pyrad = str(read_shotgun2rad_params[1])
        bwa = str(read_shotgun2rad_params[2])
        bwa_threads = int(read_shotgun2rad_params[3])
        samtools = str(read_shotgun2rad_params[4])
        mapref = str(read_shotgun2rad_params[5])
        try:
            minoverlap = float(read_shotgun2rad_params[6])
        except ValueError:
            minoverlap = 0.9
        try:
            minMAPQ = int(read_shotgun2rad_params[7])
        except ValueError:
            minMAPQ = 60
        try:
            clustsize = int(read_shotgun2rad_params[8])
        except ValueError:
            clustsize = mindepth

        """ expand ./ ~ and ../ designators in location names """

        def expander(namepath):
            if "~" in namepath:
                namepath = namepath.replace("~", os.path.expanduser("~"))
            if "../" in namepath:
                a, b = namepath.split("../")
                namepath = os.path.abspath(os.path.join(os.path.dirname(""), '..', b))
            elif "./" in namepath:
                a, b = namepath.split("./")
                namepath = os.path.abspath("") + "/" + b
            return namepath

        if WORK == "":
            WORK = os.path.abspath("") + "/"
        else:
            WORK = expander(WORK)
        if WORK[-1] != "/":
            WORK = WORK + "/"

        k = tuple('1234567')
        if options.steps:
            k = tuple(str(options.steps))

        """ shotgun2rad step 1 """
        if '1' in k:
            if not cmd_exists(pyrad):
                print "\tcannot find pyRAD, edit path in shotgun2rad_params.txt file"
                sys.exit()

            sys.stderr.write("\n\n\tshotgun2rad step 1: since WGS data is assumed to be demultiplexed, \n\t" \
                             " the analysis starts with quality filtering using pyRAD's step 2\n\t")

            pyrad_cmd = pyrad + " -p " + params_pyrad + " -s 2"
            subprocess.call(pyrad_cmd, shell=True)

        """ shotgun2rad step 2 """
        if '2' in k:
            if not cmd_exists(bwa):
                print "\tcannot find bwa, edith in shotgun2rad_params.txt"
                sys.exit()
            if not cmd_exists(samtools):
                print "\tcannot find samtools, edith in shotgun2rad_params.txt"
                sys.exit()
            if not os.path.exists(WORK + 'edits/'):
                print "\terror: could not find edits/ folder in working directory"
                sys.exit()
            if len(glob.glob(mapref + "*")) == 0:
                print "\tbwa-indexed reference not found, edith in shotgun2rad_params.txt"
            if not os.path.exists(WORK + "maps"):
                os.makedirs(WORK + "maps")
            outfolder = WORK + "maps"

            edits = glob.glob(WORK + "edits/*.edit")
            edits_sortedbysize = []

            if bwa_threads > parallel:
                parallel_s2 = 1
                bwa_threads = parallel
            else:
                parallel_s2 = parallel // bwa_threads

            sys.stderr.write("\n\n\tshotgun2rad step 2: mapping of " + `len(edits)` + " samples to " + `mapref` + "\n\t" + \
                             " Running " + `parallel_s2` + " parallel jobs; bwa-mem running with " + str(
                bwa_threads) + " threads.\n\t" + \
                             " If needed, adjust to avoid CPU and MEM limits\n\t")

            " if not only 1 sample "
            if len(edits) > 1:
                for f in edits:
                    " append files to list if not already clustered or empty"
                    if not os.path.exists(outfolder + "/" + f.split("/")[-1].replace(".edit", ".sam")):
                        size = os.stat(f)
                        if size.st_size > 0:
                            edits_sortedbysize.append(f)
                        else:
                            print "excluding " + str(f) + " file is empty"
                    else:
                        print "\tmaps/" + f.split("/")[-1].replace(".edit", ".sam") + " already exists, skipping..."

                " arranges files by decreasing size for fast SAM parsing order"
                for i in range(len(edits_sortedbysize)):
                    statinfo = os.stat(edits_sortedbysize[i])
                    edits_sortedbysize[i] = edits_sortedbysize[i], statinfo.st_size
                edits_sortedbysize.sort(key=operator.itemgetter(1), reverse=True)
                edits_sortedbysize = [i[0] for i in edits_sortedbysize]

            elif len(edits) == 1:
                f = edits
                size = os.stat(f[0])
                if size.st_size > 0:
                    edits_sortedbysize = f
                else:
                    print "excluding " + f[0] + " file is empty"
            else:
                print "\tNo .edit files found in edits/ dir."

            work_queue = multiprocessing.Queue()
            result_queue = multiprocessing.Queue()

            submitted = {}
            fileno = 1
            for edit in edits_sortedbysize:
                work_queue.put([edit, bwa, bwa_threads, mapref, samtools])
                submitted[edit] = 1
                fileno += 1

            jobs = []
            for i in range(min(submitted, parallel_s2)):
                worker = Worker(work_queue, result_queue, mapeditedreads)
                jobs.append(worker)
                worker.start()
            for j in jobs:
                j.join()


        """ shotgun2rad step 3 """
        if '3' in k:
            if not cmd_exists(muscle):
                print "\tcannot find muscle, edit path in pyRAD's param file"
                sys.exit()
            if not os.path.exists(WORK + 'maps/'):
                print "\terror: could not find maps/ folder in working directory"
                sys.exit()
            if not os.path.exists(WORK + "clust" + wclust):
                os.makedirs(WORK + "clust" + wclust)
            outfolder = WORK + "clust" + wclust

            sams = glob.glob(WORK + "maps/*.sam")
            sams_sortedbysize = []

            sys.stderr.write("\n\n\tshotgun2rad step 3: extracting loci alignments for " + `len(sams)` + " samples\n\t" + \
                             " to create .clustS.gz files. Running " + `parallel` + " parallel jobs.\n\t")

            " if not only 1 sample "
            if len(sams) > 1:
                for f in sams:
                    " append files to list if not already clustered or empty"
                    if not os.path.exists(outfolder + "/" + f.split("/")[-1].replace(".sam", ".clustS.gz")):
                        size = os.stat(f)
                        if size.st_size > 0:
                            sams_sortedbysize.append(f)
                        else:
                            print "excluding " + str(f) + " file is empty"
                    else:
                        print "\tclust" + wclust + "/" + f.split("/")[-1].replace(".sam",
                                                                                  ".clustS.gz") + " already exists, skipping..."
                " arranges files by decreasing size for fast parsing order"
                for i in range(len(sams_sortedbysize)):
                    statinfo = os.stat(sams_sortedbysize[i])
                    sams_sortedbysize[i] = sams_sortedbysize[i], statinfo.st_size
                sams_sortedbysize.sort(key=operator.itemgetter(1), reverse=True)
                sams_sortedbysize = [i[0] for i in sams_sortedbysize]

            elif len(sams) == 1:
                f = sams
                size = os.stat(f[0])
                if size.st_size > 0:
                    sams_sortedbysize = f
                else:
                    print "excluding " + f[0] + " file is empty"
            else:
                print "\tNo .sam files found in maps/ dir."

            work_queue = multiprocessing.Queue()
            result_queue = multiprocessing.Queue()

            submitted = {}
            fileno = 1
            for sam in sams_sortedbysize:
                work_queue.put([sam, outfolder, minMAPQ, minoverlap, muscle, clustsize])
                submitted[sam] = 1
                fileno += 1

            jobs = []
            for i in range(min(submitted, parallel)):
                worker = Worker(work_queue, result_queue, makeclustS)
                jobs.append(worker)
                worker.start()
            for j in jobs:
                j.join()


        """ shotgun2rad step 4 """
        if '4' in k:
            if not cmd_exists(pyrad):
                print "\tcannot find pyRAD, edit path in shotgun2rad_params.txt file"
                sys.exit()

            if not os.path.exists(WORK + "stats"):
                os.makedirs(WORK + "stats")

            pyrad_cmd = pyrad + " -p " + params_pyrad + " -s 4"
            subprocess.call(pyrad_cmd, shell=True)

        """ shotgun2rad step 5 """
        if '5' in k:
            if not cmd_exists(pyrad):
                print "\tcannot find pyRAD, edit path in shotgun2rad_params.txt file"
                sys.exit()

            if not os.path.exists(WORK + "stats"):
                os.makedirs(WORK + "stats")

            pyrad_cmd = pyrad + " -p " + params_pyrad + " -s 5"
            subprocess.call(pyrad_cmd, shell=True)

        """ shotgun2rad step 6 """
        if '6' in k:
            if not cmd_exists(pyrad):
                print "\tcannot find pyRAD, edit path in shotgun2rad_params.txt file"
                sys.exit()

            if not os.path.exists(WORK + "stats"):
                os.makedirs(WORK + "stats")

            pyrad_cmd = pyrad + " -p " + params_pyrad + " -s 6"
            subprocess.call(pyrad_cmd, shell=True)

        """ shotgun2rad step 7 """
        if '7' in k:
            if not cmd_exists(pyrad):
                print "\tcannot find pyRAD, edit path in shotgun2rad_params.txt file"
                sys.exit()

            if not os.path.exists(WORK + "stats"):
                os.makedirs(WORK + "stats")

            pyrad_cmd = pyrad + " -p " + params_pyrad + " -s 7"
            subprocess.call(pyrad_cmd, shell=True)

    """ create new shotgun2rad_params.txt file """
    if options.new_shotgun2radparams_file:
        if os.path.exists("./shotgun2rad_params.txt"):
            print "\tfile shotgun2rad_params.txt already exists"
            sys.exit()
        else:
            create_shotgun2rad_params(parser.version.split(" ")[1])


if __name__ == "__main__":
    main()

