#from __future__ import division
import random
import argparse
import sys
# bp length of mouse transcriptome in OMA: 37.914.531
# bp length of CANVA genome 2.5Mpb

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", nargs='+', default=None, help="input FASTQ filename")
parser.add_argument("-out", "--output", help="output FASTQ filename")
parser.add_argument("-f", "--fraction", type=float, help="fraction of reads to sample")
parser.add_argument("-n", "--number", type=int, help="number of reads to sample")
parser.add_argument("-s", "--sample", type=int, help="number of output files to write", default=1)
args = parser.parse_args()

if args.fraction and args.number:
   sys.exit("give either a fraction or a number, not both")

if not args.fraction and not args.number:
   sys.exit("you must give either a fraction or a number")

print("counting records....")
with open(args.input[0]) as input:
    num_lines = sum([1 for line in input])
total_records = int(num_lines / 4)

if args.fraction:
    args.number = int(total_records * args.fraction)

print("sampling " + str(args.number) + " out of " + str(total_records) + " records")

output_sequence_sets = []
output_file_left = []
if len(args.input) > 1:
    output_file_right = []
for i in range(args.sample):
    output_sequence_sets.append(set(random.sample(range(total_records + 1), args.number)))
    #output_file = args.input[0].split("/")[-1].split(".")[0]
    output_file = args.output
    output_file_left.append(open(output_file + "_0_" + str(i) + ".fq", "w"))
    if len(args.input) > 1:
        output_file_right.append(open(output_file + "_1_" + str(i) + ".fq", "w"))

initial_length = 0
sampling_length = 0

record_number = 0
with open(args.input[0]) as read_input:
    for line1 in read_input:
        line2 = read_input.readline()
        initial_length += len(line2)
        line3 = read_input.readline()
        line4 = read_input.readline()
        for i, output in enumerate(output_file_left):
            if record_number in output_sequence_sets[i]:
                    output.write(line1)
                    output.write(line2)
                    output.write(line3)
                    output.write(line4)
                    sampling_length += len(line2)
        record_number += 1

if len(args.input) > 1:
    record_number = 0
    with open(args.input[1]) as read_input:
        for line1 in read_input:
            line2 = read_input.readline()
            line3 = read_input.readline()
            line4 = read_input.readline()
            for i, output in enumerate(output_file_right):
                if record_number in output_sequence_sets[i]:
                        output.write(line1)
                        output.write(line2)
                        output.write(line3)
                        output.write(line4)
                        sampling_length += len(line2)
            record_number += 1


#output[0].close()
# if len(args.input) > 1:
#     output[1].close()
print("The mean length of all reads is {} and the mean length of the subsampled reads is {}".format(initial_length/total_records, sampling_length/args.number))
print("The sum length of all reads is {} and the sum length of the subsampled reads is {}".format(initial_length, sampling_length))
print("done!")
