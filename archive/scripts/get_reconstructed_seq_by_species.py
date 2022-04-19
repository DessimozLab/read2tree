import sys
import os
import getopt
import glob
import pandas as pd

from Bio import SeqIO
from tables import *
from Bio.SeqIO.FastaIO import FastaWriter


def read_seq_records(folder):
    out_dic = {}
    for file in glob.glob(os.path.join(folder, "*.fa")):
        sp_name = os.path.basename(file).split("_")[0]
        out_dic[sp_name] = {rec.id: rec for rec in list(SeqIO.parse(file, "fasta"))}
    return out_dic

def read_sc_file(file):
    tmp = pd.read_csv(file)
    return [t['gene_id']+"_"+t['og']+"_"+t['og'] for i,t in tmp.iterrows()]


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:s:h", ["mapping_folder=", "sc_file="])
    except getopt.GetoptError as e:
        print(str(e))
        print('get_reconstructed_seq_by_species.py -m <mapping_folder> -s <sc_file>')
        sys.exit(2)

    mapping_folder = None
    sc_file = None

    for opt, arg in opts:
        if opt == '-h':
            print('get_reconstructed_seq_by_species.py -m <mapping_folder> -s <sc_file>')
            sys.exit()
        elif opt in ("-m", "--mapping_folder"):
            mapping_folder = arg
            if mapping_folder[-1] is not "/":
                mapping_folder += "/"
        elif opt in ("-s", "--sc_file"):
            sc_file = arg
        else:
            assert False, "unhandled option"

    all_records = read_seq_records(mapping_folder)
    selected_seq = [all_records[idx[0:5]][idx] for idx in read_sc_file(sc_file)]
    print(selected_seq)
    file_name = mapping_folder.split("03_mapping_")[-1].split("/")[0]+"_consensus.fa"
    handleF = open(file_name, "w")
    writer = FastaWriter(handleF, wrap=None)
    writer.write_file(selected_seq)
    handleF.close()

if __name__ == "__main__":
    main()
