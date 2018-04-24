#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=4:00:0
#$ -pe smp 4
#$ -l tmpfs=150G
#$ -j y
#$ -N FALSP
#$ -wd /home/ucbpdvd/Scratch/output
srr=SRR5270425
folder=FALSP
cd $TMPDIR
echo 'In '$TMPDIR
~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr:0:6}/$srr/$srr.sra ./
echo 'Finished Download'
parallel-fastq-dump -s *.sra -t 4 -O . --split-files
echo 'Finished Split Files'
mkdir /home/ucbpdvd/Scratch/avian/reads/$folder
mv *_1.fastq /home/ucbpdvd/Scratch/avian/reads/$folder/$folder\_1.fq
mv *_2.fastq /home/ucbpdvd/Scratch/avian/reads/$folder/$folder\_2.fq
echo 'Finished moving files'