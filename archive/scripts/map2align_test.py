import sys
import os
import getopt
import glob

from Bio import SeqIO
from zoo.wrappers.treebuilders import Fasttree
from tables import *
from Bio import AlignIO
from zoo.wrappers.aligners import Mafft
from Bio.SeqIO import FastaIO


from zoo.seq_utils.utils import concatenate


def get_coverage(og):
    return (len(og[-1].seq)-og[-1].seq.count('X'))/len(og[-1].seq)

def perform_mapping(DIR_MAPPING, FILE_OGS):
    og_dict = {}
    '''read in og with aa seq'''
    og = list(SeqIO.parse(FILE_OGS, "fasta"))
    for record in og:
        key = record.description.split(" | ")[-1]
        if key in og_dict:
            ids = [rec.id for rec in og_dict[key]]
            if record.id not in ids:
                og_dict[key].append(record)
        else:
            og_dict[key] = []
            og_dict[key].append(record)


    # parse the mapped reads to ogs to dictionary
    all_dict = {}
    for file in glob.glob(DIR_MAPPING + "*.fa"):
        og_name = file.split("_")[-1].split(".")[0]
        og = og_dict[og_name]

        # change ids to species names
        for i, record in enumerate(og):
            s = record.id[0:5]
            record.id = s
        all_dict[og_name] = og

    OG_OUT = DIR_MAPPING + 'origin_og/'
    if not os.path.exists(OG_OUT):
        os.makedirs(OG_OUT)

    for key, item in all_dict.items():
        file_name = OG_OUT + key + ".fa"
        fasta_out = FastaIO.FastaWriter(open(file_name, "w"), wrap=None)
        fasta_out.write_file(item)

    print("FINISHED PARSING OGs!")
    return all_dict

def read_alignments(folder):
    align_list = []
    for filename in glob.glob(folder+"*.phy"):
        # input_handle = open(filename, "rU")
        align_list.append(AlignIO.read(filename, "phylip-relaxed"))
    print("FINISHED READING ALIGNMENTS!")
    return align_list

def perform_alignment(all_dict, DIR_MAPPING):
    align_dict = {}
    align_list = []
    counter = 0
    for key, value in all_dict.items():
        mafft_wrapper = Mafft(value, datatype="PROTEIN")
        mafft_wrapper.options.options['--localpair'].set_value(True)
        mafft_wrapper.options.options['--maxiterate'].set_value(1000)
        alignment = mafft_wrapper()
        align_dict[key] = alignment
        align_list.append(alignment)
        counter += 1
        if counter % 50 == 0:
            print('{} of {} alignments done'.format(counter, len(all_dict)))

    ALIGN_OUT = DIR_MAPPING + 'origin_align/'

    if not os.path.exists(ALIGN_OUT):
        os.makedirs(ALIGN_OUT)
    print("WRITING ALIGNMENT FILES INTO: {}!".format(ALIGN_OUT))
    for key, value in align_dict.items():
        output_handle = open(ALIGN_OUT + key + ".phy", "w")
        AlignIO.write(value, output_handle, "phylip")
    print("FINISHED ALIGNMENTS!")
    return align_list

def concatenate_alignment(align_list, DIR_MAPPING):
    ALIGN_OUT = DIR_MAPPING + 'origin_align/'
    concat_align = concatenate(align_list)

    output_handle = open(ALIGN_OUT + "CONCAT.phy", "w")
    AlignIO.write(concat_align, output_handle, "phylip")
    print("FINISHED CONCATINATION!")
    return concat_align

def build_tree(concat_align, DIR_MAPPING):

    fasttree_wrapper = Fasttree(concat_align, datatype="PROTEIN")
    tree = fasttree_wrapper()
    print("FINISHED TREE INFERENCE!")
    with open(DIR_MAPPING+"original_tree.nwk", "w") as text_file:
        text_file.write("{};".format(tree))
    print("Resulting tree: {}".format(tree))
    return tree

def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:o:a:t:h", ["mapping_folder=", "ortholog_file=", "alignmnet_folder="])
    except getopt.GetoptError as e:
        print(str(e))
        print('map2align_test.py -m <mapping_folder> -o <ortholog_file> -a <alignmnet_folder>')
        sys.exit(2)

    mapping_folder = None
    ortholog_file = None
    alignment_folder = None

    for opt, arg in opts:
        if opt == '-h':
            print('map2align_test.py -m <mapping_folder> -o <ortholog_file> -a <alignmnet_folder>')
            sys.exit()
        elif opt in ("-m", "--mapping_folder"):
            mapping_folder = arg
            if mapping_folder[-1] is not "/":
                mapping_folder += "/"
        elif opt in ("-a", "--alignmnet_folder"):
            alignment_folder = arg
            if alignment_folder[-1] is not "/":
                alignment_folder += "/"
        elif opt in ("-o", "--ortholog_folder"):
            ortholog_file = arg
        else:
            assert False, "unhandled option"


    if alignment_folder:
        align = read_alignments(alignment_folder)
    else:
        mapping = perform_mapping(mapping_folder, ortholog_file)
        align = perform_alignment(mapping, mapping_folder)

    concatenation = concatenate_alignment(align, mapping_folder)
    build_tree(concatenation, mapping_folder)


if __name__ == "__main__":
    main()
