import sys, re
import os
import getopt
import glob
import pyopa


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from zoo.wrappers.treebuilders import Fasttree
from tables import *
from Bio import AlignIO
from zoo.wrappers.aligners import Mafft
from Bio.SeqIO import FastaIO


from zoo.seq_utils.utils import concatenate

defaults = pyopa.load_default_environments()
envs = defaults['environments']
env = envs[515]


def find_best_translation_by_similarity(mapped_sequences, reference_og_sequence, exclude_species=["not_me"]):
    '''Given a list of sequences that are derived from mapped reads to multiple seq of a OG
    we find the best corresponding mapped seq by comparing it with a representative sequence of the original OG using
    pyopa local alignment and return the sequence with its highest score!'''
    best_score = 0
    best_sequence = None
    s1 = pyopa.Sequence(str(reference_og_sequence.seq))
    for record in mapped_sequences:
        if record.id[0:5] not in exclude_species:
            # print(record)
            frames = [record.seq[i:].translate(table='Standard', stop_symbol='*', to_stop=True, cds=False, gap="N") for
                      i in range(3)]
            best_seq_idx = 0
            for i, seq in enumerate(frames):
                s2 = pyopa.Sequence(str(seq))
                # calculating local and global scores for the given sequences
                local_double = pyopa.align_double(s1, s2, env)
                # print('Local score: %f' % local_double[0])
                if local_double[0] > best_score:
                    best_score = local_double[0]
                    best_seq_idx = i
                    best_sequence = SeqRecord(frames[best_seq_idx], id="simul", description=record.description,
                                              name=record.name)
                    # print(best_sequence)

    return best_sequence


def get_coverage(og):
    return (len(og[-1].seq)-og[-1].seq.count('X'))/len(og[-1].seq)


def perform_mapping(mapping_dir, og_files, threshold=1, exclude_species=["none"]):
    og_dict = {}
    '''read in og with aa seq'''
    og = list(SeqIO.parse(og_files, "fasta"))
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
    for file in glob.glob(mapping_dir + "*.fa"):
        og_name = file.split("_")[-1].split(".")[0]
        og = og_dict[og_name]

        # change ids to species names
        for i, record in enumerate(og):
            s = record.id[0:5]
            record.id = s

        # find the best representative seq based by mapping
        mapping = list(SeqIO.parse(file, "fasta"))
        best_translated_seq = find_best_translation_by_similarity(mapping, og[0], exclude_species=exclude_species)
        if best_translated_seq is not None:
            og.append(best_translated_seq)
            if get_coverage(og) <= threshold:
                all_dict[og_name] = og

    if threshold is not 1:
        OG_OUT = mapping_dir + 'og' + str(threshold) + "/"
    elif exclude_species[0] is not "none":
        OG_OUT = mapping_dir + 'og_without' + "_".join(exclude_species) + "/"
    else:
        OG_OUT = mapping_dir + 'og/'

    if not os.path.exists(OG_OUT):
        os.makedirs(OG_OUT)

    for key, item in all_dict.items():
        file_name = OG_OUT + key + ".fa"
        fasta_out = FastaIO.FastaWriter(open(file_name, "w"), wrap=None)
        fasta_out.write_file(item)

    print("FINISHED OG RECONSTRUCTION!")
    return all_dict


def read_alignments(folder):
    align_list = []
    for filename in glob.glob(folder+"*.phy"):
        # input_handle = open(filename, "rU")
        align_list.append(AlignIO.read(filename, "phylip-relaxed"))
    print("FINISHED READING ALIGNMENTS!")
    return align_list


def perform_alignment(all_dict, mapping_dir, threshold=1, exclude_species=["none"]):
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

    if threshold is not 1:
        ALIGN_OUT = mapping_dir + 'align' + str(threshold) + "/"
    elif exclude_species[0] is not "none":
        ALIGN_OUT = mapping_dir + 'align_without' + "_".join(exclude_species) + "/"
    else:
        ALIGN_OUT = mapping_dir + 'align/'

    if not os.path.exists(ALIGN_OUT):
        os.makedirs(ALIGN_OUT)
    print("WRITING ALIGNMENT FILES INTO: {}!".format(ALIGN_OUT))

    for key, value in align_dict.items():
        output_handle = open(ALIGN_OUT + key + ".phy", "w")
        AlignIO.write(value, output_handle, "phylip")
    print("FINISHED ALIGNMENTS!")

    return align_list


def concatenate_alignment(align_list, mapping_dir, threshold=1, exclude_species=["none"]):

    if threshold is not 1:
        ALIGN_OUT = mapping_dir + 'align' + str(threshold) + "/"
    elif exclude_species[0] is not "none":
        ALIGN_OUT = mapping_dir + 'align_without' + "_".join(exclude_species) + "/"
    else:
        ALIGN_OUT = mapping_dir + 'align/'

    concat_align = concatenate(align_list)

    output_handle = open(ALIGN_OUT + "CONCAT.phy", "w")
    AlignIO.write(concat_align, output_handle, "phylip")
    print("FINISHED CONCATINATION!")
    return concat_align


def build_tree(concat_align, mapping_dir, threshold=1, exclude_species=["none"]):

    if threshold is not 1:
        tree_out = mapping_dir + 'tree' + str(threshold) + ".nwk"
    elif exclude_species[0] is not "none":
        tree_out = mapping_dir + 'tree_without' + "_".join(exclude_species) + ".nwk"
    else:
        tree_out = mapping_dir + 'tree.nwk'

    fasttree_wrapper = Fasttree(concat_align, datatype="PROTEIN")
    tree = fasttree_wrapper()
    print("FINISHED TREE INFERENCE!")
        with open(tree_out, "w") as text_file:
            text_file.write("{};".format(tree))
    print("Resulting tree: {}".format(tree))
    return tree


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:o:a:t:e:h", ["mapping_folder=", "ortholog_file=", "alignmnet_folder=", "threshold=", "exclude_species="])
    except getopt.GetoptError as e:
        print(str(e))
        print('map2align.py -m <mapping_folder> -o <ortholog_file> -a <alignmnet_folder> -t <threshold> -e <exclude_species>')
        sys.exit(2)

    mapping_folder = None
    ortholog_file = None
    alignment_folder = None
    threshold = 1
    exclude_species = ["none"]

    for opt, arg in opts:
        if opt == '-h':
            print('map2align.py -m <mapping_folder> -o <ortholog_file> -a <alignmnet_folder> -t <threshold> -e <exclude_species>')
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
        elif opt in ("-t", "--threshold"):
            threshold = float(arg)
        elif opt in ("-e", "--exclude_species"):
            exclude_species = arg.split(",")
        else:
            assert False, "unhandled option"


    if alignment_folder:
        align = read_alignments(alignment_folder)
    else:
        mapping = perform_mapping(mapping_folder, ortholog_file, threshold=threshold, exclude_species=exclude_species)
        align = perform_alignment(mapping, mapping_folder, threshold=threshold, exclude_species=exclude_species)

    concatenation = concatenate_alignment(align, mapping_folder, threshold=threshold, exclude_species=exclude_species)
    build_tree(concatenation, mapping_folder, threshold=threshold, exclude_species=exclude_species)


if __name__ == "__main__":
    main()
