#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=0:10:0
#$ -pe smp 1
#$ -j y
#$ -N rm_FALSP
#$ -wd /home/ucbpdvd/Scratch/output
rm -r /home/ucbpdvd/Scratch/avian/reads/FALSP