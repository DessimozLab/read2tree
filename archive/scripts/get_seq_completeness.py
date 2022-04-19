import sys
import os
import getopt
import glob

from Bio import SeqIO
from tables import *
from read2tree.stats.SeqCompleteness import SeqCompleteness


def read_seq_records(file):
    return list(SeqIO.parse(file, "fasta"))


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:r:h", ["mapping_folder=", "reference_folder="])
    except getopt.GetoptError as e:
        print(str(e))
        print('get_seq_completeness.py -m <mapping_folder> -r <reference_folder>')
        sys.exit(2)

    mapping_folder = None
    reference_folder = None

    for opt, arg in opts:
        if opt == '-h':
            print('get_seq_completeness.py -m <mapping_folder> -r <reference_folder>')
            sys.exit()
        elif opt in ("-m", "--mapping_folder"):
            mapping_folder = arg
            if mapping_folder[-1] is not "/":
                mapping_folder += "/"
        elif opt in ("-r", "--reference_folder"):
            reference_folder = arg
            if reference_folder[-1] is not "/":
                reference_folder += "/"
        else:
            assert False, "unhandled option"
    if reference_folder and mapping_folder:
        for file in glob.glob(reference_folder+"/*.fa"):
            species = os.path.basename(file).split("_")[0]
            print(species)
            ref_records = read_seq_records(file)
            mapping_file = os.path.join(mapping_folder,os.path.basename(file).split(".")[0]+"_consensus.fa")
            if os.path.exists(mapping_file):
                map_records = read_seq_records(mapping_file)
                seqC = SeqCompleteness(ref_records)
                seqC.get_seq_completeness(map_records)
                seqC.write_seq_completeness(os.path.join(mapping_folder, species + "_OGs_sc.txt"))

if __name__ == "__main__":
    main()
