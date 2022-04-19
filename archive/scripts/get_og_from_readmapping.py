import sys
import getopt
from Bio import SeqIO
from Bio.SeqIO import FastaIO

def get_ogs(mapped_reads, og_data):
    og_data_names = [record.name for record in og_data]
    og_data_ogs = [record.description.split("| ")[-1] for record in og_data]
    list_of_ogs = {}
    for record in mapped_reads:
        if record.name in og_data_names:
            og_index = og_data_names.index(record.name)
            og_name = og_data[og_index].description.split("| ")[-1]
            indices = [i for i, x in enumerate(og_data_ogs) if x == og_name]
            seq_to_write = [og_data[i] for i in indices]
            record.seq = record.seq.upper()
            record.id = "SRR400661_" + record.id
            record.name = "SRR400661_" + record.name
            record.description = "SRR400661_" + record.description
            seq_to_write.append(record)
            list_of_ogs[og_name] = seq_to_write
    return list_of_ogs



def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "r:d:o:", ["mapped_reads=", "ref_data=", "out_folder="])
    except getopt.GetoptError as e:
        print(str(e))
        print('concat_alignments.py -r <mapped_reads> -d <ref_data> -o <out_folder>')
        sys.exit(2)

    mapped_reads = None
    ref_data = None
    out_folder = None

    for opt, arg in opts:
        if opt == '-h':
            print('concat_alignments.py -r <mapped_reads> -d <ref_data> -o <out_folder>')
            sys.exit()
        elif opt in ("-r", "--reads"):
            mapped_reads = arg
        elif opt in ("-d", "--ref_data"):
            ref_data = arg
        elif opt in ("-o", "--out_folder"):
            out_folder = arg
        else:
            assert False, "unhandled option"

    read_mappings = list(SeqIO.parse(mapped_reads, "fasta"))
    og_data = list(SeqIO.parse(ref_data, "fasta"))

    if out_folder[-1] is not "/":
        out_folder += "/"

    list_of_ogs = get_ogs(read_mappings, og_data)
    if list_of_ogs is not None:
        for og in list_of_ogs:
            file_name = out_folder + og + ".fasta"
            fasta_out = FastaIO.FastaWriter(open(file_name, "w"), wrap=None)
            fasta_out.write_file(list_of_ogs[og])


if __name__ == "__main__":
    main()