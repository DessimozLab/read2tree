import sys
import os
import getopt
import glob
from Bio import AlignIO, SeqIO

from zoo.seq_utils.utils import concatenate

<<<<<<< Updated upstream
def concatenate_alignments(folder):
=======
def concatenate_alignments(folder, min_taxa=0):
>>>>>>> Stashed changes
    all_og_alignments = []
    all_og_align_pos = {}
    start = 0
    for f in glob.glob(folder+'*.phy'):
        used_ogs = 0
        if os.path.getsize(f) > 0:
            try:
                msa = AlignIO.read(f, "phylip-relaxed")
            except ValueError:
                msa = AlignIO.read(f, "fasta")
        #for record in msa:
        #    record.id = record.id[0:5]
        #msa[-1].id = "CANAL"
            if len(msa) >= min_taxa:
                print(f)
                used_ogs =+ 1
                all_og_alignments.append(msa)
            #all_og_align_pos[f] = [start, start + len(record.seq)]
            #start = len(record.seq) + 1
    con_alignment = concatenate(all_og_alignments)
    print('OGs used: {}!'.format(used_ogs))
    return con_alignment


def main():

    try:
      opts, args = getopt.getopt(sys.argv[1:], "f:m:o:", ["folder=", "min_taxa=", "out_file="])
    except getopt.GetoptError as e:
        print(str(e))
        print('concat_alignments.py -f <folder> -m <min_taxa> -o <out_file>')
        sys.exit(2)

    seq_folder = None
    out_file = None
    min_taxa = 0

    for opt, arg in opts:
        if opt == '-h':
            print('concat_alignments.py -f <folder> -m <min_taxa> -o <out_folder> -d')
            sys.exit()
        elif opt in ("-f", "--folder"):
            seq_folder = arg
        elif opt in ("-o", "--out_file"):
            out_file = arg
        elif opt in ("-m", "--min_taxa"):
            min_taxa = int(arg)
        else:
            assert False, "unhandled option"



    if seq_folder[-1] is not "/":
        seq_folder += "/"
    
    if min_taxa > 0:
      out_file = out_file+"_"+str(min_taxa)+".phy"

    alignment = concatenate_alignments(seq_folder, min_taxa=min_taxa)
    if alignment is not None:
        align_output = open(out_file, "w")
        AlignIO.write(alignment, align_output, "phylip-relaxed")
        align_output.close()

if __name__ == "__main__":
    main()
