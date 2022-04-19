#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/avian/lsf_out/down_GLYSP.o%J
#BSUB -e /scratch/beegfs/weekly/ddylus/avian/lsf_out/down_GLYSP.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J down_GLYSP
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -M 2000000
srr=SRR3115005
speciesid=GLYSP
module add Utility/aspera_connect/3.7.4.147727
module add UHTS/Analysis/sratoolkit/2.8.2.1
source activate r2t
mkdir /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
ascp -v -QT -k1 -l100M -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr:0:6}/$srr/$srr.sra ./
echo 'Finished download'
fastq-dump --split-files --gzip $srr.sra
echo 'Finished getting fastq from sra and split files'
mv *\_1.* $speciesid\_1.fq.gz
mv *\_2.* $speciesid\_2.fq.gz
echo 'Finished moving files'