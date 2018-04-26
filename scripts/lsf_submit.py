import os
import sys
import getopt
import glob
import time
import subprocess
import pandas as pd
import numpy as np

def get_sra_dic(df):
    all_index = []
    for organism in list(set(df.Organism)):
        if organism is not np.nan:
            data_by_organism = df[df.Organism == organism]
            print(len(data_by_organism))
            if len(data_by_organism) == 1:  # there is only one datapoint for this organism lets collect it anyway
                if data_by_organism.MBases.values[
                    0] > 100:  # for sure select the ones that have more than 1X coverage assuming a 1GB genome
                    all_index.append(data_by_organism.index[0])
                    print(data_by_organism['MBases'])
            else:
                # Check if there is transcriptomic data available
                if 'TRANSCRIPTOMIC' in data_by_organism.LibrarySource:
                    data_o_t = data_by_organism[df.LibrarySource == 'TRANSCRIPTOMIC']
                    x = data_by_organism.iloc[(data_by_organims['MBases'] - 1000).abs().argsort()[:1]]
                    if x.MBases.values[
                        0] > 1000:  # for sure select the ones that have more than 10X coverage assuming a 1GB genome
                        all_index.append(x.index[0])
                elif 'GENOMIC' in data_by_organism.LibrarySource:
                    x = data_by_organism.iloc[(data_by_organims['MBases'] - 10000).abs().argsort()[:1]]
                    if x.MBases.values[
                        0] > 10000:  # for sure select the ones that have more than 10X coverage assuming a 1GB genome
                        all_index.append(x.index[0])
    #                 else: #  run closest subset sum algorithm to select multiple ones that we are able to merge
    #                     data_by_organism_with_most_technology =

    sra_dic = {}
    for index, row in x.iterrows():
        sra_dic[row.Organism] = row.Run
    return sra_dic


def get_download_string(species_id, sra, se_pe='PAIRED'):
    if 'ERR' in sra and se_pe is 'PAIRED':
        download = """#!/bin/bash
#BSUB -o {}.o%J
#BSUB -e {}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J {}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
srr={}
speciesid={}
module add Utility/aspera_connect/3.7.4.147727
source activate r2t
mkdir /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
ascp -v -QT -k1 -l100M -P33001 -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${srr:0:6}/002/$srr/$srr\_1.fastq.gz .
ascp -v -QT -k1 -l100M -P33001 -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${srr:0:6}/002/$srr/$srr\_2.fastq.gz .
echo 'Finished download'
gunzip -c $srr\_1.fastq.gz > $speciesid\_1.fq
gunzip -c $srr\_2.fastq.gz > $speciesid\_2.fq
echo 'Finished moving files'""".format(species_id, species_id, species_id, sra, species_id)
    if 'ERR' in sra and se_pe is 'SINGLE':
        download = """#!/bin/bash
#BSUB -o {}.o%J
#BSUB -e {}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J {}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
srr={}
speciesid={}
module add Utility/aspera_connect/3.7.4.147727
source activate r2t
mkdir /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
ascp -v -QT -k1 -l100M -P33001 -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${srr:0:6}/002/$srr/$srr.fastq.gz .
echo 'Finished download'
gunzip -c $srr.fastq.gz > $speciesid\_1.fq
echo 'Finished moving files'""".format(species_id, species_id, species_id, sra, species_id)
    elif 'SRR' in sra and se_pe is 'PAIRED':
        download = """#!/bin/bash
#BSUB -o {}.o%J
#BSUB -e {}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J {}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
srr={}
speciesid={}
module add Utility/aspera_connect/3.7.4.147727
source activate r2t
mkdir /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
ascp -v -QT -k1 -l100M -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr:0:6}/$srr/$srr.sra ./
echo 'Finished download'
parallel-fastq-dump -s *.sra -t 4 -O . --split-files
echo 'Finished getting fastq from sra and split files'
mv *\_1.fastq $speciesid\_1.fq
mv *\_2.fastq $speciesid\_2.fq
echo 'Finished moving files'""".format(species_id, species_id, species_id, sra, species_id)
    elif 'SRR' in sra and se_pe is 'SINGLE':
        download = """#!/bin/bash
#BSUB -o {}.o%J
#BSUB -e {}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J {}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
srr={}
speciesid={}
module add Utility/aspera_connect/3.7.4.147727
source activate r2t
mkdir /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
ascp -v -QT -k1 -l100M -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr:0:6}/$srr/$srr.sra ./
echo 'Finished download'
parallel-fastq-dump -s *.sra -t 4 -O .
echo 'Finished getting fastq from sra'
mv *.fastq $speciesid\_1.fq
echo 'Finished moving files'""".format(species_id, species_id, species_id, sra, species_id)

    text_file = open('down_py_script.sh', "w")
    text_file.write(download)
    text_file.close()

    return 'down_py_script.sh'


def get_r2t_string(species_id, reference, se_pe='PAIRED', read_type='short'):
    if se_pe is 'PAIRED' and read_type is 'short':
        job_string = """#!/bin/bash
#BSUB -o {}.o%J
#BSUB -e {}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J {}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
source activate r2t
reads=/scratch/beegfs/weekly/ddylus/avian/reads/{}
cd /scratch/beegfs/weekly/ddylus/avian/r2t/
source activate r2t
python -W ignore ~/opt/read2tree/bin/read2tree --standalone_path /scratch/beegfs/weekly/ddylus/avian/marker_genes/ --dna_reference /scratch/beegfs/weekly/ddylus/avian/eukaryotes.cdna.fa --reads $reads/{}_1.fq $reads/{}_2.fq --output_path /scratch/beegfs/weekly/ddylus/avian/r2t/ --single_mapping {} --threads 4 --min_species 8""".format(species_id, species_id, species_id, species_id, reference)
    elif se_pe is 'SINGLE' and read_type is 'short':
        job_string = """#!/bin/bash
#BSUB -o {}.o%J
#BSUB -e {}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J {}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
source activate r2t
reads=/scratch/beegfs/weekly/ddylus/avian/reads/{}
cd /scratch/beegfs/weekly/ddylus/avian/r2t/
source activate r2t
python -W ignore ~/opt/read2tree/bin/read2tree --standalone_path /scratch/beegfs/weekly/ddylus/avian/marker_genes/ --dna_reference /scratch/beegfs/weekly/ddylus/avian/eukaryotes.cdna.fa --reads $reads/{}_1.fq --output_path /scratch/beegfs/weekly/ddylus/avian/r2t/ --single_mapping {} --threads 4 --min_species 8""".format(
        species_id, species_id, species_id, reference)
    elif se_pe is 'SINGLE' and read_type is 'long':
        job_string = """#!/bin/bash
#BSUB -o {}.o%J
#BSUB -e {}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J {}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000
source activate r2t
reads=/scratch/beegfs/weekly/ddylus/avian/reads/{}
cd /scratch/beegfs/weekly/ddylus/avian/r2t/
source activate r2t
python -W ignore ~/opt/read2tree/bin/read2tree --standalone_path /scratch/beegfs/weekly/ddylus/avian/marker_genes/ --dna_reference /scratch/beegfs/weekly/ddylus/avian/eukaryotes.cdna.fa --reads $reads/{}_1.fq --output_path /scratch/beegfs/weekly/ddylus/avian/r2t/ --single_mapping {} --threads 4 --min_species 8 --read_type long""".format(
            species_id, species_id, species_id, reference)

    text_file = open('r2t_py_script.sh', "w")
    text_file.write(job_string)
    text_file.close()

    return 'r2t_py_script.sh'


def get_rm_string(species_id):
    rm = """#!/bin/bash
#BSUB -o {}.o%J
#BSUB -e {}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J {}
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -M 1000000
#$ -wd /home/ucbpdvd/Scratch/output
rm -r /scratch/beegfs/weekly/ddylus/avian/reads/{}""".format(species_id, species_id, species_id, species_id)

    text_file = open('rm_py_script.sh', "w")
    text_file.write(rm)
    text_file.close()
    return 'rm_py_script.sh'


def get_five_letter_species_id(species):
    tmp = species.split(" ")
    new_id = tmp[0][:3].upper() + tmp[1][:2].upper()
    return new_id


def is_species_mapped(species_id, output):
    if os.path.exists(os.path.join(output, '03_mapping_'+species_id+'_1')):
        mapped_speciesid = os.path.join(output, '03_mapping_'+species_id+'_1')
        files = [f for f in glob.glob(os.path.join(mapped_speciesid, "*.fa"))]
        if files:
            if len(files) == 10:
                return True
            else:
                return False
        else:
            return False
    else:
        return False

def output_shell(line):
    """
    Save output of shell line that has pipes
    taken from: https://stackoverflow.com/questions/7389662/link-several-popen-commands-with-pipes
    :param line:
        :return:
"""
    try:
        shell_command = subprocess.Popen(line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    except OSError:
        return None
    except ValueError:
        return None

    (output, err) = shell_command.communicate()
    shell_command.wait()
    if shell_command.returncode != 0:
        print("Shell command failed to execute")
        print(line)
        return None

    return output

def run_lsf(sra_dic, output_speciesid):
    num_job_cycles = 0
    rm_job_id_idx = 0
    rm_job_id = []
    for species, sra in sra_dic.items():
        species_id = get_five_letter_species_id(species)
        if not is_species_mapped(species_id, output_speciesid):  # check whether the mapping already exists
            print('Submitting species {} with species id {}!'.format(species, species_id))
            if num_job_cycles < 3:  # only run three jobs, then submit the jobs with dependency that files are again deleted
                # Set up download string
                job_string = get_download_string(species_id, sra[0], se_pe=sra[1])

                # Open a pipe to the bsub command.
                p_download = output_shell('bsub ' + job_string)
                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = p_download.decode("utf-8").split(" ")[2]

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_speciesid, '02_ref_dna/*.fa')):
                    # Set up r2t string
                    r2t_job_string = get_r2t_string(species_id, ref, se_pe=sra[1], read_type=sra[2])

                    # Open a pipe to the bsub command.
                    #output_r2t, input_r2t = Popen('bsub -hold_jid {}'.format(jobid))
                    p_download = output_shell('bsub -hold_jid {} {}'.format(jobid, r2t_job_string))

                    # Append jobid of r2t
                    r2t_jobids.append(p_download.decode("utf-8").split(" ")[2])

                    time.sleep(0.1)

                # Set up r2t string
                rm_job_string = get_rm_string(species_id)

                # Open a pipe to the bsub command.
                #output_rm, input_rm = Popen('bsub -hold_jid {}'.format(','.join(r2t_jobids)))
                p_rm = output_shell('bsub -hold_jid {} {}'.format(','.join(r2t_jobids), rm_job_string))

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(p_rm.decode("utf-8").split(" ")[2])
                time.sleep(0.1)
            else:  # this part should ensure that on the scratch never more than 3 downloads are available
                # Set up download string
                job_string = get_download_string(species_id, sra[0], se_pe=sra[1])

                # Open a pipe to the bsub command.
                p_download = output_shell('bsub -hold_jid {} {}'.format(rm_job_id[rm_job_id_idx], job_string))

                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = p_download.decode("utf-8").split(" ")[2]

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_speciesid, '02_ref_dna/*.fa')):
                    # Set up r2t string
                    r2t_job_string = get_r2t_string(species_id, ref, se_pe=sra[1], read_type=sra[2])

                    # Open a pipe to the bsub command.
                    # output_r2t, input_r2t = Popen('bsub -hold_jid {}'.format(jobid))
                    p_download = output_shell('bsub -hold_jid {} {}'.format(jobid, r2t_job_string))

                    # Append jobid of r2t
                    r2t_jobids.append(p_download.decode("utf-8").split(" ")[2])

                    time.sleep(0.1)

                # Set up r2t string
                rm_job_string = get_rm_string(species_id)

                # Open a pipe to the bsub command.
                # output_rm, input_rm = Popen('bsub -hold_jid {}'.format(','.join(r2t_jobids)))
                p_rm = output_shell('bsub -hold_jid {} {}'.format(','.join(r2t_jobids), rm_job_string))

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(p_rm.decode("utf-8").split(" ")[2])
                rm_job_id_idx += 1
                time.sleep(0.1)
            num_job_cycles += 1

def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:o:", ["sra_file=", "out_speciesid="])
    except getopt.GetoptError as e:
        print(str(e))
        print('sge_submit.py -s <sra_file> -o <out_speciesid>')
        sys.exit(2)

    sra_file = None
    out_speciesid = None


    for opt, arg in opts:
        if opt == '-h':
            print('sge_submit.py -s <sra_file> -m <min_taxa> -o <out_speciesid> -d')
            sys.exit()
        elif opt in ("-s", "--sra_file"):
            sra_file = arg
        elif opt in ("-o", "--out_speciesid"):
            out_speciesid = arg
        else:
            assert False, "unhandled option"

    # df = pd.read_csv(sra_file, sep='\t')
    sra_dic = {'Zosterops minutus': ['SRR2145255', 'SINGLE', 'short'],
 'Caloperdix oculeus': ['SRR2952863', 'PAIRED', 'short'],
 'Modulatrix stictigula': ['SRR2968792', 'PAIRED', 'short'],
 'Francolinus bicalcaratus': ['SRR2959742', 'PAIRED', 'short']}

    run_lsf(sra_dic, out_speciesid)

if __name__ == "__main__":
    main()
