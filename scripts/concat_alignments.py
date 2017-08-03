import sys
import getopt
import glob
from Bio import AlignIO, SeqIO

from zoo.seq_utils.utils import concatenate

def concatenate_alignments(folder):
        all_og_alignments = []
    all_og_align_pos = {}
    start = 0
    for f in glob.glob(folder+'*.phy'):
        msa = AlignIO.read(f, "phylip-relaxed")
        for record in msa:
            record.id = record.id[0:5]
        msa[-1].id = "CANAL"
        all_og_alignments.append(msa)
        all_og_align_pos[f] = [start, start + len(record.seq)]
        start = len(record.seq) + 1
    con_alignment = concatenate(all_og_alignments)
    return con_alignment


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:o:", ["folder=", "out_file="])
    except getopt.GetoptError as e:
        print(str(e))
        print 'concat_alignments.py -f <folder> -o <out_file>'
        sys.exit(2)

    seq_folder = None
    out_file = None


    for opt, arg in opts:
        if opt == '-h':
            print 'concat_alignments.py -f <folder> -m <min_taxa> -o <out_folder> -d'
            sys.exit()
        elif opt in ("-f", "--folder"):
            seq_folder = arg
        elif opt in ("-o", "--out_file"):
            out_file = arg
        else:
            assert False, "unhandled option"



    if seq_folder[-1] is not "/":
        seq_folder += "/"


    alignment = concatenate_alignments(seq_folder)
    if alignment is not None:
        align_output = open(out_file, "w")
        AlignIO.write(alignment, align_output, "phylip-relaxed")
        align_output.close()



if __name__ == "__main__":
    main()