from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SeqUtils import GC
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys
import getopt


def gc(seq1):
    seq = Seq(seq1)
    return GC(seq)


def transcribe(seq1):
    seq = Seq(seq1)
    return seq.transcribe()


def reverse_complement(seq1):
    seq = Seq(seq1)
    return seq.reverse_complement()


def calc_nbases(seq):
    nBases = seq.count('n') + seq.count('N')
    return nBases


def is_valid(seq, seqType):
    dna = "ACGT"
    rna = "ACGU"
    protein = "ABCDEFGHIKLMNPQRSTVWXYZ"
    if seqType.upper() == "PROTEIN":
        for letter in seq:
            if letter.upper() not in protein:
                return False
        else:
            return True
    if seqType.upper() == "DNA":
        for letter in seq:
            if letter.upper() not in dna:
                return False
        else:
            return True
    if seqType.upper() == "RNA":
        for letter in seq:
            if letter.upper() not in rna:
                return False
        else:
            return True


def filter_nBases(seq):
    newseq = []
    for x in seq:
        if not (x == 'n' or x == "N"):
            newseq.append(x)
    return ''.join(newseq)


def seq_alignment(seq1, seq2, o=""):
    alignment = pairwise2.align.globalxx(seq1, seq2)
    if o != "":
        alignmentFile = open(o, 'w')
        alignmentFile.write(str(alignment))
        alignmentFile.close()
        return "Alignment file successfully created"
    else:
        return alignment


def seq_alignment_files(file1, file2, output=""):
    try:
        seqFile1 = list(SeqIO.parse(file1, "fasta"))[0]
        seqFile2 = list(SeqIO.parse(file2, "fasta"))[0]
    except IOError:
        print("Files not exist")
        return
    seq1 = seqFile1.seq
    seq2 = seqFile2.seq
    result = seq_alignment(seq1, seq2, output)
    return result


def online_alignment(seq, output=""):
    result_handle = NCBIWWW.qblast("blastn", "nt", seq)
    blast_record = NCBIXML.read(result_handle)
    threshold = 0.01
    if output != "":
        outputFile = open(output, 'w')
        for x in blast_record.alignments:
            outputFile.write(str(x.title))
            outputFile.write(str(x.length))
            for y in x.hsps:
                if y.expect < threshold:
                    outputFile.write(str(y.score))
                    outputFile.write(str(y.bits))
                    outputFile.write(str(y.expect))
                    outputFile.write(str(y.num_alignments))
                    outputFile.write(str(y.identities))
                    outputFile.write(str(y.positives))
                    outputFile.write(str(y.gaps))
                    outputFile.write(str(y.strand))
                    outputFile.write(str(y.frame))
                    outputFile.write(str(y.query))
                    outputFile.write(str(y.query_start))
                    outputFile.write(str(y.match))
                    outputFile.write(str(y.sbjct))
                    outputFile.write(str(y.sbjct_start))
            for x in blast_record.descriptions:
                outputFile.write(str(x.title))
                outputFile.write(str(x.score))
                outputFile.write(str(x.e))
                outputFile.write(str(x.num_alignments))
        outputFile.write(str(blast_record.multiple_alignment))
        outputFile.close()
        return "Alignment wrote on file successfully"
    else:
        for x in blast_record.alignments:
            print(x.title)
            print(x.length)
            for y in x.hsps:
                if y.expect < threshold:
                    print(y.score)
                    print(y.bits)
                    print(y.expect)
                    print(y.num_alignments)
                    print(y.identities)
                    print(y.positives)
                    print(y.gaps)
                    print(y.strand)
                    print(y.frame)
                    print(y.query)
                    print(y.query_start)
                    print(y.match)
                    print(y.sbjct)
                    print(y.sbjct_start)
        for x in blast_record.descriptions:
            print(x.title)
            print(x.score)
            print(x.e)
            print(x.num_alignments)
        print(blast_record.multiple_alignment)
        return "End of Blast fetching"


def merge_fasta(out, *files):
    f = ""
    file = ""
    if out != "":
        for i in files:
            try:
                f = open(i, 'r')
            except IOError:
                print("file not exist")

            f2 = open(out, 'a')
            for j in f:
                f2.write(j)
        return "Written successfully"
    elif out == "":
        for x in files:
            try:
                file = SeqIO.parse(x, "fasta")
            except IOError:
                print(x, "file not exist")

            for line in file:
                print(line)
        return "Merging Done"


# def merge_fasta_write(out, *files):
#     for i in files:
#         try:
#             f = open(i, 'r')
#         except IOError:
#             print("file not exist")
#             return
#         f2 = open(out, 'a')
#         for j in f:
#             for i in j:
#                 f2.write(i)
#
#     return "Merging Saved Successfully to the file"


def convert_to_fasta(file):
    fileFasta = file.split(".gbk")
    # print(fileFasta)
    faa = f"{fileFasta[0]}.fasta"
    SeqIO.convert(file, "genbank", faa, "fasta")
    return "File converted successfully"


def getOpt():
    opts, args = getopt.getopt(sys.argv[1:], "o:", "help")

    i = 0
    argsDict = {}
    li = ["gc", "calc_nbases", "filter_nBases", "transcribe", "reverse_complement", "convert_to_fasta", "is_valid",
          "seq_alignment", "seq_alignment_files", "online_alignment", "merge_fasta"]
    # print(len(args), args[0])
    if len(args) > 1 and args[0] in li:
        while i < len(args):

            if args[i] in ["gc", "calc_nbases", "filter_nBases", "transcribe", "reverse_complement",
                           "convert_to_fasta"]:
                value = []
                value.append(args[i + 1])
                argsDict[args[i]] = value
                i += 2
            elif args[i] in ["is_valid"]:
                isValid = args[i + 1: i + 3]
                argsDict[args[i]] = isValid
                i += 3
            elif args[i] in ["seq_alignment", "seq_alignment_files"]:
                seq = args[i + 1: i + 3]
                argsDict[args[i]] = seq
                # i = i + 3
                if len(opts) == 1:
                    lis = args[i + 1:i + 3]
                    lis.append(opts[0][1])
                    argsDict[args[i]] = lis
                    i += 4
                else:
                    i += 3

            elif args[i] in ["online_alignment"]:
                lis = []
                seq = args[i + 1]
                argsDict[args[i]] = seq
                if len(opts) == 1:
                    lis.append(args[i + 1])
                    lis.append(opts[0][1])
                    argsDict[args[i]] = lis
                    i += 3
                else:
                    i += 2
            elif args[i] in ["merge_fasta"]:
                argsDict[args[i]] = []
                fa = []
                counter = 0
                for x in args:
                    if x not in li:
                        fa.append(x)
                        counter += 1
                argsDict[args[i]] = fa
                i += counter

    elif len(args) == 1 and args[0] in li:
        print("Parameters missing")
        sys.exit(0)
    elif args[0] not in li:
        print("Incorrect Command")
        sys.exit(0)
    # print(args, opts)
    for k, v in argsDict.items():
        # print(k, v)
        if k in ["help"]:
            sys.exit()
        if k in ["gc", "calc_nbases", "filter_nBases", "transcribe", "reverse_complement", "convert_to_fasta"]:
            if len(v) <= 0:
                print("Parameters number not satisfied")
                sys.exit(0)
            else:
                if k == "gc":
                    if len(v) == 1:
                        print(gc(v[0]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
                elif k == "calc_nbases":
                    if len(v) == 1:
                        print(calc_nbases(v[0]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
                elif k == "filter_nBases":
                    if len(v) == 1:
                        print(filter_nBases(v[0]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
                elif k == "transcribe":
                    if len(v) == 1:
                        print(transcribe(v[0]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
                elif k == "reverse_complement":
                    if len(v) == 1:
                        print(reverse_complement(v[0]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
                elif k == "convert_to_fasta":
                    if len(v) == 1:
                        print(convert_to_fasta(v[0]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
        elif k in ["is_valid", "seq_alignment", "seq_alignment_files", "online_alignment", "merge_fasta"]:
            if len(v) <= 1:
                print("Parameters number not satisfied")
                sys.exit(0)
            else:
                if k == "is_valid":
                    if len(v) == 2:
                        print(is_valid(v[0], v[1]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)

                elif k == "seq_alignment":
                    if len(v) == 2 and len(opts) == 0:
                        print(seq_alignment(v[0], v[1]))
                    elif len(v) == 3 and len(opts) == 1:
                        print(seq_alignment(v[0], v[1], opts[0][1]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)

                elif k == "seq_alignment_files":
                    if len(v) == 2 and len(opts) == 0:
                        print(seq_alignment_files(v[0], v[1]))
                    elif len(v) == 3 and len(opts) == 1:
                        print(seq_alignment_files(v[0], v[1], opts[0][1]))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
                elif k == "online_alignment":
                    # print(v, "bla", v[0], len(v), len(argsDict.values()), argsDict.values(),len(opts))
                    # print(len(argsDict.values()))
                    if len(v) == 2 and len(opts) == 1:
                        # print("with file")
                        print(online_alignment(v[0], opts[0][1]))
                    elif len(argsDict.values()) == 1 and len(opts) == 0:
                        # print("without")
                        print(online_alignment(v))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
                elif k == "merge_fasta":
                    mergeV = v[1:]
                    if len(opts) == 1 and len(mergeV) >= 2:
                        print(merge_fasta(opts[0][1], *mergeV))
                    elif len(mergeV) >= 2 and len(opts) == 0:
                        print(merge_fasta("", *mergeV))
                    else:
                        print("Parameters number not satisfied")
                        sys.exit(0)
                else:
                    print("Parameters number not satisfied")
                    sys.exit(0)
        else:
            print("Unexpected Error")
            sys.exit(0)


getOpt()
# x = merge_fasta("outt.txt", "myseq.fasta", "ls_orchid.fasta")
# print(x)
