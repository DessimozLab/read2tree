import os
import sys
import getopt
import glob
import time
import subprocess
import pandas as pd


def get_name_to_id(df):
    name_to_id = {}
    for i,r in df.iterrows():
        name_to_id[r['Organism']] = r['id5letter']
    return name_to_id

def get_sra_dic(df, name_to_id):
    sra_dic = {}
    for sp, idx in name_to_id.items():
        subset = df[df.Organism == sp].sort_values(by=['MBases'], ascending=False)
        sra_dic[sp] = [idx]
        sra_dic[sp].append([r['Run'] for i,r in subset.iterrows()])
    return sra_dic


def get_download_string_ena(species_id, sra, se_pe='PAIRED'):
    print(sra)
    sra_string = ''
    for i in sra:
        sra_string += '\"'+i+'\"'
        sra_string += ' '
    download = """#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=10:00:0
#$ -pe smp 1
#$ -l tmpfs=100G
#$ -j y
#$ -N down_%s
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
speciesid=%s
source activate r2t
mkdir /home/ucbpdvd/Scratch/avian/reads/$speciesid
reads=/home/ucbpdvd/Scratch/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /home/ucbpdvd/Scratch/avian/reads/$speciesid
declare -a sra_all=(%s)
for sra in "${sra_all[@]}"
do
    echo $sra
    ~/.aspera/connect/bin/ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${sra:0:6}/00${sra: -1}/$sra/$sra\_1.fastq.gz .
    ~/.aspera/connect/bin/ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${sra:0:6}/00${sra: -1}/$sra/$sra\_2.fastq.gz .
    echo 'Finished $sra'
    if [ ! -s $sra\_1.fastq.gz ] && [ ! -s $sra\_2.fastq.gz ]
    then
        echo $sra
        ~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/${sra:0:3}/${sra:0:6}/$sra/$sra.sra .
        fastq-dump --split-files --gzip $sra.sra
    fi
done

find . -name "*\_1.*" | sort -V | xargs cat > $speciesid\_1.fq.gz
find . -name "*\_2.*" | sort -V | xargs cat > $speciesid\_2.fq.gz
python -W ignore /home/ucbpdvd/opt/read2tree/scripts/sample_reads.py --coverage 5 --genome_len 1000000000 --reads $speciesid\_1.fq.gz $speciesid\_2.fq.gz

for sra in "${sra_all[@]}"
do
    rm $sra*
done
echo 'Finished moving files'""" % (species_id, species_id, sra_string.rstrip())
    text_file = open('down_py_script.sh', "w")
    text_file.write(download)
    text_file.close()

    return 'down_py_script.sh'


def get_r2t_string(species_id, reference, se_pe='PAIRED', read_type='short'):
    if se_pe is 'PAIRED' and read_type is 'short':
        job_string = """#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=16:00:0
#$ -pe smp 4
#$ -l tmpfs=140G
#$ -j y
#$ -N r2t_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
reads=/home/ucbpdvd/Scratch/avian/reads/{species_id}
cd /home/ucbpdvd/Scratch/avian/r2t/
source activate r2t
if [ -s $reads/{species_id}\_1.fa ] && [ -s $reads/{species_id}\_2.fa ]
then
    python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fa $reads/{species_id}_2.fa \
--output_path /home/ucbpdvd/Scratch/avian/r2t/ --single_mapping {reference} \
--threads 4 --min_species 0 --read_type short
elif [ -s $reads/{species_id}\_1.fq.gz ]  && [ -s $reads/{species_id}\_2.fq.gz ]
then 
	python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fq.gz $reads/{species_id}_2.fq.gz \
--output_path /home/ucbpdvd/Scratch/avian/r2t/ --single_mapping {reference} \
--threads 4 --min_species 0 --read_type short --check_mate_pairing
fi""".format(species_id=species_id,reference=reference)

    text_file = open('r2t_py_script.sh', "w")
    text_file.write(job_string)
    text_file.close()

    return 'r2t_py_script.sh'


def get_rm_string(species_id):
    rm = """#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=0:10:0
#$ -pe smp 1
#$ -j y
#$ -N rm_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
rm -r \
/home/ucbpdvd/Scratch/avian/reads/{species_id}""".format(species_id=species_id)

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


def run_sge(sra_dic, output_folder):
    num_job_cycles = 0
    rm_job_id_idx = 0
    rm_job_id = []
    for species, sra in sra_dic.items():
        species_id = sra[0]
        sra_ids = sra[-1]
        # check whether the mapping already exists
        mapped_species = get_species_mapped(species_id, output_folder)
        number_species = len([f for f in glob.glob(os.path.join(output_folder, '02_ref_dna/*.fa'))])
        if len(mapped_species) < number_species:
            print('Submitting species {} with species id {} and mapping {}!'.format(species, species_id,
                                                                                    len(mapped_species)))
            if num_job_cycles < 3:  # only run three jobs, then submit the jobs with dependency that files are again deleted
                # Set up download string
                job_string = get_download_string_ena(species_id, sra_ids,
                                                 se_pe='PAIRED')

                # Open a pipe to the qsub command.
                p_download = output_shell('qsub ' + job_string)
                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = p_download.decode("utf-8").split(" ")[2]

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_folder,
                                                  '02_ref_dna/*.fa')):
                    if os.path.basename(ref) not in mapped_species:
                        # Set up r2t string
                        r2t_job_string = get_r2t_string(
                            species_id, ref, se_pe='PAIRED', read_type='short')

                        # Open a pipe to the qsub command.
                        # output_r2t, input_r2t = Popen('qsub -hold_jid {}'.format(jobid))
                        p_download = output_shell('qsub -hold_jid {} {}'
                                                  .format(jobid, r2t_job_string))

                        # Append jobid of r2t
                        r2t_jobids.append(p_download.decode("utf-8").split(" ")[2])

                        time.sleep(0.1)

                # Set up r2t string
                rm_job_string = get_rm_string(species_id)

                # Open a pipe to the qsub command.
                # output_rm, input_rm = Popen('qsub -hold_jid {}'.format(','.join(r2t_jobids)))
                p_rm = output_shell(
                    'qsub -hold_jid {} {}'
                    .format(','.join(r2t_jobids), rm_job_string))

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(p_rm.decode("utf-8").split(" ")[2])
                time.sleep(0.1)
            else:  # this part should ensure that on the scratch never more than 3 downloads are available
                # Set up download string
                job_string = get_download_string_ena(species_id, sra_ids,
                                                 se_pe='PAIRED')

                # Open a pipe to the qsub command.
                p_download = output_shell(
                    'qsub -hold_jid {} {}'.format(rm_job_id[rm_job_id_idx],
                                                  job_string))

                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = p_download.decode("utf-8").split(" ")[2]

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_folder,
                                                  '02_ref_dna/*.fa')):
                    if os.path.basename(ref) not in mapped_species:
                        # Set up r2t string
                        r2t_job_string = get_r2t_string(
                            species_id, ref, se_pe='PAIRED', read_type='short')

                        # Open a pipe to the qsub command.
                        # output_r2t, input_r2t = Popen('qsub -hold_jid {}'.format(jobid))
                        p_download = output_shell('qsub -hold_jid {} {}'
                                                  .format(jobid, r2t_job_string))

                        # Append jobid of r2t
                        r2t_jobids.append(p_download.decode("utf-8").split(" ")[2])

                        time.sleep(0.1)

                # Set up r2t string
                rm_job_string = get_rm_string(species_id)

                # Open a pipe to the qsub command.
                # output_rm, input_rm = Popen('qsub -hold_jid {}'.format(','.join(r2t_jobids)))
                p_rm = output_shell(
                    'qsub -hold_jid {} {}'.format(','.join(r2t_jobids),
                                                  rm_job_string))

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(p_rm.decode("utf-8").split(" ")[2])
                rm_job_id_idx += 1
                time.sleep(0.1)
            num_job_cycles += 1


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Submit jobs given an sra file""")
    parser.add_argument('--sra_file', default=None, required=True,
                        help='[Default is none] SRA-file with ')
    parser.add_argument('--output_folder', default=None, required=True,
                            help='[Default is none] Folder of r2t run.')

    conf = parser.parse_args()
    df = pd.read_csv(conf.sra_file)
    df_illumina_paired = df[(df.Platform == 'ILLUMINA')]
    print('Make sure to set the folders of reads, '
          'cluster ouput and location of run correctly!')

    print('We are selecting {} / {} illumina paired entries.'
          .format(len(set(df_illumina_paired.Organism)), len(df)))

    # df = pd.read_csv(sra_file, sep='\t')
    name_to_id = get_name_to_id(df_illumina_paired)
    sra_dic = get_sra_dic(df_illumina_paired, name_to_id)

    run_sge(sra_dic, conf.output_folder)
