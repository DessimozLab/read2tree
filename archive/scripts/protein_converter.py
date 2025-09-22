import sys
from Bio import SeqIO

# Get input and output file paths from command-line arguments
# Daniel Paiva Agustinho
input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, "r") as input_handle:
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            protein_seq = record.seq.translate()

            # Extract the entire original header
            original_header = record.description

            # Create a new sequence record with the original header
            protein_seq = SeqIO.SeqRecord(
                protein_seq, id=record.id, description=original_header
            )
            SeqIO.write(protein_seq, output_handle, "fasta")
print("done",str(output_file))
