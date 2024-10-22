

import read2tree
from read2tree.main import main
from read2tree._utils import exe_name

import sys

print("start run_r2t 2  223  3 ")
main(sys.argv[1:], exe_name=exe_name(), desc="descr")


# --step 1marker  --standalone_path marker_genes  --dna_reference dna_ref.fa --output_path output  --debug
# --step 2map --standalone_path marker_genes  --dna_reference dna_ref.fa --reads /work/FAC/FBM/DBC/cdessim2/read2tree/v2_test/t1/reads_20/ERR7350657__2.fastq.gz  --output_path output --debug --threads 1

print("finish run_r2t   ")



a=1