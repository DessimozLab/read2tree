import os
import re
import sys
import getopt
import glob
import time
import subprocess
import pandas as pd


def get_download_string(species_id, sra, layout='PAIRED'):
    sra_string = ''
    for i in sra:
        sra_string += '\"'+i+'\"'
        sra_string += ' '
    download = """#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/metazoan/lsf/down_%s.o%sJ
#BSUB -e /scratch/beegfs/weekly/ddylus/metazoan/lsf/down_%s.e%sJ
#BSUB -u david.dylus@unil.ch
#BSUB -J down_%s
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -M 2000000
speciesid=%s
module add Bioinformatics/Software/vital-it
module add Utility/aspera_connect/3.7.4.147727
module add UHTS/Analysis/sratoolkit/2.8.2.1
reads=/scratch/beegfs/weekly/ddylus/metazoan/reads/$speciesid
mkdir /scratch/beegfs/weekly/ddylus/metazoan/reads/$speciesid
echo Created read $speciesid
cd /scratch/beegfs/weekly/ddylus/metazoan/reads/$speciesid
declare -a sra_all=(%s)
if [ "%s" == "PAIRED" ]
then
    for sra in "${sra_all[@]}"
    do
        ascp -QT -l 300m -P33001 -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${sra:0:6}/00${sra: -1}/$sra/$sra\_1.fastq.gz .
        ascp -QT -l 300m -P33001 -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${sra:0:6}/00${sra: -1}/$sra/$sra\_2.fastq.gz .
        echo 'Finished download'
        # rm $sra.sra
        if [ ! -s $sra\_1.fastq.gz ] && [ ! -s $sra\_1.fastq.gz ]
        then
            echo ----- USING NCBI DOWNLOAD -----
            ascp -v -QT -k1 -l100M -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh  anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/${sra:0:3}/${sra:0:6}/$sra/$sra.sra .
            fastq-dump --split-files --gzip $sra.sra
        fi
    done
    
    find . -name "*\_1.*" | sort -V | xargs cat > $speciesid\_1.fq.gz
    find . -name "*\_2.*" | sort -V | xargs cat > $speciesid\_2.fq.gz
    python -W ignore /scratch/beegfs/weekly/ddylus/opt/read2tree/scripts/sample_reads.py --coverage 20 --genome_len 50000000 --reads $speciesid\_1.fq.gz $speciesid\_2.fq.gz
    for sra in "${sra_all[@]}"
    do
        rm $sra*
    done
else
    for sra in "${sra_all[@]}"
    do
        ascp -QT -l 300m -P33001 -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${sra:0:6}/00${sra: -1}/$sra/$sra.fastq.gz .
        echo 'Finished download'
        # rm $sra.sra
        if [ ! -s "$sra.fastq.gz" ]
            then
            echo $sra
            ascp -v -QT -k1 -l100M -i /software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh  anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/${sra:0:3}/${sra:0:6}/$sra/$sra.sra .
            fastq-dump --gzip $sra.sra
        fi
    done
    find . -name "*.gz" | sort -V | xargs cat > $speciesid\_1.fq.gz
    for sra in "${sra_all[@]}"
    do
        rm $sra*
    done
fi

echo 'Finished moving files'""" % (species_id, '%', species_id, '%', species_id, species_id, sra_string.rstrip(), layout)


    text_file = open('down_py_script.sh', "w")
    text_file.write(download)
    text_file.close()

    return 'down_py_script.sh'


def get_r2t_string(species_id, reference, se_pe='PAIRED', read_type='short'):
    if se_pe in 'PAIRED' and read_type is 'short':
        job_string = """#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/metazoan/lsf/r2t_{species_id}.o%J
#BSUB -e /scratch/beegfs/weekly/ddylus/metazoan/lsf/r2t_{species_id}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J r2t_{species_id}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=12000]"
#BSUB -M 12000000
conda activate r2t
reads=/scratch/beegfs/weekly/ddylus/metazoan/reads/{species_id}
cd /scratch/beegfs/weekly/ddylus/metazoan/r2t/
if [ -s $reads/{species_id}\_1.fa ] && [ -s $reads/{species_id}\_2.fa ]
then
    python -W ignore /scratch/beegfs/weekly/ddylus/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fa $reads/{species_id}_2.fa \
--output_path . --single_mapping 02_ref_dna/{reference}
--threads 4 --min_species 0 --read_type short 
elif [ -s $reads/{species_id}\_1.fq.gz ]  && [ -s $reads/{species_id}\_2.fq.gz ]
then
    python -W ignore /scratch/beegfs/weekly/ddylus/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fq.gz $reads/{species_id}_2.fq.gz \
--output_path . --single_mapping 02_ref_dna/{reference}
--threads 4 --min_species 0 --read_type short 
fi""".format(species_id=species_id, reference=os.path.basename(reference))
    elif se_pe in 'SINGLE' and read_type is 'short':
        job_string = """#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/metazoan/lsf/r2t_{species_id}.o%J
#BSUB -e /scratch/beegfs/weekly/ddylus/metazoan/lsf/r2t_{species_id}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J r2t_{species_id}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=12000]"
#BSUB -M 12000000
conda activate r2t
reads=/scratch/beegfs/weekly/ddylus/metazoan/reads/{species_id}
cd /scratch/beegfs/weekly/ddylus/metazoan/r2t/
python -W ignore /scratch/beegfs/weekly/ddylus/opt/read2tree/bin/read2tree  \
--reads $reads/{species_id}_1.fa \
--output_path . --single_mapping 02_ref_dna/{reference} \
--threads 4 --min_species 0 --read_type short""".format(species_id=species_id, reference=os.path.basename(reference))
    elif se_pe is 'SINGLE' and read_type is 'long':
        job_string = """#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/metazoan/lsf/r2t_{species_id}.o%J
#BSUB -e /scratch/beegfs/weekly/ddylus/metazoan/lsf/r2t_{species_id}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J r2t_{species_id}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=12000]"
#BSUB -M 12000000
conda activate r2t
reads=/scratch/beegfs/weekly/ddylus/metazoan/reads/{species_id}
cd /scratch/beegfs/weekly/ddylus/metazoan/r2t/
python -W ignore /scratch/beegfs/weekly/ddylus/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fa \
--output_path . \
--single_mapping 02_ref_dna/{reference} --threads 4 --min_species 8 --read_type short \
--split_reads""".format(species_id=species_id, reference=os.path.basename(reference))
    text_file = open('r2t_py_script.sh', "w")
    text_file.write(job_string)
    text_file.close()

    return 'r2t_py_script.sh'


def get_rm_string(species_id):
    rm = """#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/metazoan/lsf/rm_%s.o%sJ
#BSUB -e /scratch/beegfs/weekly/ddylus/metazoan/lsf/rm_%s.e%sJ
#BSUB -u david.dylus@unil.ch
#BSUB -J rm_%s
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -M 1000000
rm -r /scratch/beegfs/weekly/ddylus/metazoan/reads/%s""" % (species_id, '%', species_id, '%', species_id, species_id)

    text_file = open('rm_py_script.sh', "w")
    text_file.write(rm)
    text_file.close()
    return 'rm_py_script.sh'


def get_five_letter_species_id(species):
    tmp = species.split(" ")
    new_id = tmp[0][:3].upper() + tmp[1][:2].upper()
    return new_id


def get_species_mapped(species_id, output):
    if os.path.exists(os.path.join(output, '04_mapping_'+species_id+'_1')):
        mapped_speciesid = os.path.join(output, '04_mapping_'+species_id+'_1')
        files = [os.path.basename(f).replace('_consensus.fa','.fa') for f in glob.glob(os.path.join(mapped_speciesid, "*.fa"))]
        if files:
            if len(files):
                return files
            else:
                return []
        else:
            return []
    else:
        return []


def output_shell(line):
    """
    Save output of shell line that has pipes
    taken from: https://stackoverflow.com/questions/7389662/link-several-popen-commands-with-pipes
    :param line:
        :return:
"""
    try:
        shell_command = subprocess.Popen(
            line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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


def run_lsf(sra_dic, output_r2t):
    num_job_cycles = 0
    rm_job_id_idx = 0
    rm_job_id = []
    for species, sra in sra_dic.items():
        species_id = sra[-1]
        # check whether the mapping already exists
        mapped_species = get_species_mapped(species_id, output_r2t)
        number_species = len([f for f in glob.glob(os.path.join(output_r2t, '02_ref_dna/*.fa'))])
        if len(mapped_species) < number_species:
            print('Submitting species {} with species id {} and mapping {}!'.format(species, species_id, len(mapped_species)))
            if num_job_cycles < 3:  # only run three jobs, then submit the jobs with dependency that files are again deleted
                # Set up download string
                job_string = get_download_string(species_id, sra[0], sra[1])

                # Open a pipe to the bsub command.
                p_download = output_shell('bsub < ' + job_string)
                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = re.search('<(.*)>', p_download.decode("utf-8").split(" ")[1]).group(1)

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_r2t, '02_ref_dna/*.fa')):
                    if os.path.basename(ref) not in mapped_species:
                        # Set up r2t string
                        r2t_job_string = get_r2t_string(species_id, ref, se_pe=sra[1], read_type=sra[2])

                        # Open a pipe to the bsub command.
                        #output_r2t, input_r2t = Popen('bsub -hold_jid {}'.format(jobid))
                        p_download = output_shell(
                            'bsub -w "{}" < {}'.format('done('+jobid+')', r2t_job_string))

                        # Append jobid of r2t
                        r2t_jobids.append(
                            re.search('<(.*)>', p_download.decode("utf-8").split(" ")[1]).group(1))

                        time.sleep(0.1)

                # Set up r2t string
                rm_job_string = get_rm_string(species_id)

                # Open a pipe to the bsub command.
                #output_rm, input_rm = Popen('bsub -hold_jid {}'.format(','.join(r2t_jobids)))
                r2t_jobids_done = ['done(' + r + ')' for r in r2t_jobids]
                p_rm = output_shell(
                    'bsub -w "{}" < {}'.format('&&'.join(r2t_jobids_done), rm_job_string))

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(re.search('<(.*)>', p_rm.decode("utf-8").split(" ")[1]).group(1))
                time.sleep(0.1)
            else:  # this part should ensure that on the scratch never more than 3 downloads are available
                # Set up download string
                job_string = get_download_string(species_id, sra[0], sra[1])

                # Open a pipe to the bsub command.
                p_download = output_shell(
                    'bsub -w "{}" < {}'.format('done('+rm_job_id[rm_job_id_idx]+')', job_string))

                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = re.search('<(.*)>', p_download.decode("utf-8").split(" ")[1]).group(1)

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_r2t, '02_ref_dna/*.fa')):
                    if os.path.basename(ref) not in mapped_species:
                        # Set up r2t string
                        r2t_job_string = get_r2t_string(species_id, ref, se_pe=sra[1], read_type=sra[2])

                        # Open a pipe to the bsub command.
                        # output_r2t, input_r2t = Popen('bsub -hold_jid {}'.format(jobid))
                        p_download = output_shell(
                            'bsub -w "{}" < {}'.format('done('+jobid+')', r2t_job_string))

                        # Append jobid of r2t
                        r2t_jobids.append(
                            re.search('<(.*)>', p_download.decode("utf-8").split(" ")[1]).group(1))

                        time.sleep(0.1)

                # Set up r2t string
                rm_job_string = get_rm_string(species_id)

                # Open a pipe to the bsub command.
                # output_rm, input_rm = Popen('bsub -hold_jid {}'.format(','.join(r2t_jobids)))
                r2t_jobids_done = ['done('+r+')' for r in r2t_jobids]
                p_rm = output_shell(
                    'bsub -w "{}" < {}'.format('&&'.join(r2t_jobids_done), rm_job_string))

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(re.search('<(.*)>', p_rm.decode("utf-8").split(" ")[1]).group(1))
                rm_job_id_idx += 1
                time.sleep(0.1)
            num_job_cycles += 1


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:o:", ["sra_file=", "output_r2t="])
    except getopt.GetoptError as e:
        print(str(e))
        print('sge_submit.py -s <sra_file> -o <output_r2t>')
        sys.exit(2)

    sra_file = None
    output_r2t = None

    for opt, arg in opts:
        if opt == '-h':
            print('sge_submit.py -s <sra_file> -o <output_r2t> -d')
            sys.exit()
        elif opt in ("-s", "--sra_file"):
            sra_file = arg
        elif opt in ("-o", "--output_r2t"):
            output_r2t = arg
        else:
            assert False, "unhandled option"

    df = pd.read_csv(sra_file, sep=',')
    # species = ['XENTR', 'STRPU', 'ATTCE', 'LOTGI', 'TRIAD', 'AMPQE']
    species = []
    sra_dic = {}
    for i, row in df.iterrows():
        if row['oma_ids'] not in species:
            sra_dic[row['species_ids']] = [[row['sra'].rstrip()], row['layout'].rstrip(), 'short', row['oma_ids'].rstrip()]

    run_lsf(sra_dic, output_r2t)


if __name__ == "__main__":
    main()
