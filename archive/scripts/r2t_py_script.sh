#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/avian/lsf_out/r2t_GLYSP.o%J
#BSUB -e /scratch/beegfs/weekly/ddylus/avian/lsf_out/r2t_GLYSP.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J r2t_GLYSP
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
source activate r2t
reads=/scratch/beegfs/weekly/ddylus/avian/reads/GLYSP
cd /scratch/beegfs/weekly/ddylus/avian/r2t/
python -W ignore /scratch/beegfs/monthly/ddylus/opt/read2tree/bin/read2tree --standalone_path /scratch/beegfs/weekly/ddylus/avian/marker_genes/ --dna_reference /scratch/beegfs/weekly/ddylus/avian/eukaryotes.cdna.fa --reads $reads/GLYSP_1.fq.gz $reads/GLYSP_2.fq.gz --output_path /scratch/beegfs/weekly/ddylus/avian/r2t/ --single_mapping /scratch/beegfs/weekly/ddylus/avian/r2t/02_ref_dna/MELGA_OGs.fa --threads 4 --min_species 8