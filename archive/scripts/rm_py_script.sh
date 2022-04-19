#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/avian/lsf_out/rm_GLYSP.o%J
#BSUB -e /scratch/beegfs/weekly/ddylus/avian/lsf_out/rm_GLYSP.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J rm_GLYSP
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -M 1000000
rm -r /scratch/beegfs/weekly/ddylus/avian/reads/GLYSP