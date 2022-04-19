import os
import re
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


def get_download_string(species_id, sra, layout='PAIRED'):
    sra_string = ''
    for i in sra:
        sra_string += '\"'+i+'\"'
        sra_string += ' '
    download = """#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/avian/lsf/down_%s.o%sJ
#BSUB -e /scratch/beegfs/weekly/ddylus/avian/lsf/down_%s.e%sJ
#BSUB -u david.dylus@unil.ch
#BSUB -J down_%s
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -M 2000000
speciesid=%s
module add Utility/aspera_connect/3.7.4.147727
module add UHTS/Analysis/sratoolkit/2.8.2.1
reads=/scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
mkdir /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
echo Created read $speciesid
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
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
    python -W ignore /scratch/beegfs/weekly/ddylus/opt/read2tree/scripts/sample_reads.py --coverage 5 --genome_len 1000000000 --reads $speciesid\_1.fq.gz $speciesid\_2.fq.gz
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
#BSUB -o /scratch/beegfs/weekly/ddylus/avian/lsf/r2t_{species_id}.o%J
#BSUB -e /scratch/beegfs/weekly/ddylus/avian/lsf/r2t_{species_id}.e%J
#BSUB -u david.dylus@unil.ch
#BSUB -J r2t_{species_id}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=1000]"
#BSUB -M 4000000
reads=/scratch/beegfs/weekly/ddylus/avian/reads/{species_id}
cd /scratch/beegfs/weekly/ddylus/avian/r2t/
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

    text_file = open('r2t_py_script.sh', "w")
    text_file.write(job_string)
    text_file.close()

    return 'r2t_py_script.sh'


def get_rm_string(species_id):
    rm = """#!/bin/bash
#BSUB -o /scratch/beegfs/weekly/ddylus/avian/lsf/rm_%s.o%sJ
#BSUB -e /scratch/beegfs/weekly/ddylus/avian/lsf/rm_%s.e%sJ
#BSUB -u david.dylus@unil.ch
#BSUB -J rm_%s
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -M 1000000
rm -r /scratch/beegfs/weekly/ddylus/avian/reads/%s""" % (species_id, '%', species_id, '%', species_id, species_id)

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


def run_lsf(sra_dic, output_folder, max_jobs):
    num_total_submissions = 0
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
                job_string = get_download_string(species_id, sra_ids)

                # Open a pipe to the bsub command.
                p_download = output_shell('bsub < ' + job_string)
                num_total_submissions += 1
                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = re.search('<(.*)>', p_download.decode("utf-8").split(" ")[1]).group(1)

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_folder, '02_ref_dna/*.fa')):
                    # Set up r2t string
                    if os.path.basename(ref) not in mapped_species:
                        r2t_job_string = get_r2t_string(species_id, ref, se_pe='PAIRED', read_type='short')

                        # Open a pipe to the bsub command.
                        #output_r2t, input_r2t = Popen('bsub -hold_jid {}'.format(jobid))
                        p_download = output_shell(
                            'bsub -w "{}" < {}'.format('done('+jobid+')', r2t_job_string))
                        num_total_submissions += 1

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
                    'bsub -w "{}" < {}'.format(' && '.join(r2t_jobids_done), rm_job_string))
                num_total_submissions += 1

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(re.search('<(.*)>', p_rm.decode("utf-8").split(" ")[1]).group(1))
                time.sleep(0.1)
            else:  # this part should ensure that on the scratch never more than 3 downloads are available
                # Set up download string
                job_string = get_download_string(species_id, sra_ids)

                # Open a pipe to the bsub command.
                p_download = output_shell(
                    'bsub -w "{}" < {}'.format('done('+rm_job_id[rm_job_id_idx]+')', job_string))
                num_total_submissions += 1
                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = re.search('<(.*)>', p_download.decode("utf-8").split(" ")[1]).group(1)

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_folder, '02_ref_dna/*.fa')):
                    if os.path.basename(ref) not in mapped_species:
                        # Set up r2t string
                        r2t_job_string = get_r2t_string(species_id, ref, se_pe='PAIRED', read_type='short')

                        # Open a pipe to the bsub command.
                        # output_r2t, input_r2t = Popen('bsub -hold_jid {}'.format(jobid))
                        p_download = output_shell(
                            'bsub -w "{}" < {}'.format('done('+jobid+')', r2t_job_string))
                        num_total_submissions += 1

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
                    'bsub -w "{}" < {}'.format(' && '.join(r2t_jobids_done), rm_job_string))
                num_total_submissions += 1

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(re.search('<(.*)>', p_rm.decode("utf-8").split(" ")[1]).group(1))
                rm_job_id_idx += 1
                time.sleep(0.1)
            num_job_cycles += 1
        if num_total_submissions > max_jobs:
            break


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Submit jobs given an sra file""")
    parser.add_argument('--sra_file', default=None, required=True,
                        help='[Default is none] SRA-file with ')
    parser.add_argument('--output_folder', default=None, required=True,
                            help='[Default is none] Folder of r2t run.')
    parser.add_argument('--max_jobs', type=int, default=1000,
                        help='[Default is 1000] Number of jobs.')


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

    run_lsf(sra_dic, conf.output_folder, conf.max_jobs)


