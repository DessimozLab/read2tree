#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/avian/lsf_out/down_GLYSP.o%J
#BSUB -e /scratch/beegfs/weekly/ddylus/avian/lsf_out/down_GLYSP.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J down_GLYSP
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
srr=SRR3115005
speciesid=GLYSP
module add Utility/aspera_connect/3.7.4.147727
source activate r2t
mkdir /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
ascp -v -QT -k1 -l100M -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr:0:6}/$srr/$srr.sra ./
echo 'Finished download'
parallel-fastq-dump -s *.sra -t 4 -O . --split-files --tmpdir .
echo 'Finished getting fastq from sra and split files'
mv *\_1.fastq $speciesid\_1.fq
mv *\_2.fastq $speciesid\_2.fq
echo 'Finished moving files'