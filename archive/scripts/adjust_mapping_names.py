import sys
import shutil
import os
import getopt
import glob

from Bio import SeqIO
from tables import *

def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:r:h", ["mapping_folder=", "reference_folder="])
    except getopt.GetoptError as e:
        print(str(e))
        print('get_seq_completeness.py -m <mapping_folder>')
        sys.exit(2)

    mapping_folder = None

    for opt, arg in opts:
        if opt == '-h':
            print('get_seq_completeness.py -m <mapping_folder>')
            sys.exit()
        elif opt in ("-m", "--mapping_folder"):
            mapping_folder = arg
            if mapping_folder[-1] is not "/":
                mapping_folder += "/"
        else:
            assert False, "unhandled option"

    if mapping_folder:
        for file in glob.glob(mapping_folder + "/*.fa"):
            if "_OGs" not in os.path.basename(file):
                species_name = os.path.basename(file).split("_")[0]
                new_file_name = species_name + "_OGs_consensus.fa"
                shutil.move(file, os.path.join(mapping_folder, new_file_name))
        for file in glob.glob(mapping_folder + "/*cov.txt"):
            if "_OGs" not in os.path.basename(file):
                species_name = os.path.basename(file).split("_")[0]
                new_file_name = species_name + "_OGs_cov.txt"
                shutil.move(file, os.path.join(mapping_folder, new_file_name))

if __name__ == "__main__":
    main()
