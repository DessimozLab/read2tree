#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=8:00:0
#$ -pe smp 4
#$ -l tmpfs=150G
#$ -j y
#$ -N r2t_FALSP
#$ -wd /home/ucbpdvd/Scratch/output
reads=/home/ucbpdvd/Scratch/avian/reads/FALSP
cd /home/ucbpdvd/Scratch/avian/r2t/
source activate r2t
python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree --standalone_path /home/ucbpdvd/Scratch/avian/marker_genes/ --dna_reference /home/ucbpdvd/Scratch/avian/eukaryotes.cdna.fa --reads $reads/FALSP_1.fq $reads/FALSP_2.fq --output_path /home/ucbpdvd/Scratch/avian/r2t/ --single_mapping /home/ucbpdvd/Scratch/avian/r2t/02_ref_dna/MELGA_OGs.fa --threads 4 --min_species 8