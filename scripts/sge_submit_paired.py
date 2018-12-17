
import os
import sys
import getopt
import glob
import time
import subprocess
import itertools
import numpy as np
from collections import OrderedDict
import pandas as pd


def get_name_to_id(df):
    name_to_id = {}
    index = 0
    for ogr in sorted(set(df.Organism)):
        #     print(ogr)
        if len(ogr.split(" ")) > 2:
            new_id = ogr.split(" ")[0][0:3].upper()
            + ogr.split(" ")[2][0:2].upper()
        elif len(ogr.split(" ")) > 1:
            new_id = ogr.split(" ")[0][0:3].upper()
            + ogr.split(" ")[1][0:2].upper()
        else:
            new_id = ogr.split(" ")[0][0:3].upper()+'sp'

        if new_id in name_to_id.values():
            use_id = new_id[0:4]+str(index)
    #         print(use_id)
            index += 1
        else:
            index = 0
            use_id = new_id
        name_to_id[ogr] = use_id
    return name_to_id


def get_sra_dic(df, name_to_id, start=0, end=0):
    all_index = []
    sra_dic = {}
    for organism in list(set(df.Organism)):
        if organism is not np.nan:
            data_by_organism = df[df.Organism == organism]
    #         print(len(data_by_organism))
            if len(data_by_organism) == 1:  # there is only one datapoint for this organism lets collect it anyway
                # for sure select the ones that have more than 1X coverage assuming a 1GB genome
                if data_by_organism.MBases.values[0] > 100:
                    sra_dic[organism] = [name_to_id[organism],
                                         data_by_organism.MBases.values[0]]
                    all_index.append([data_by_organism.index[0]])
                    sra_dic[organism].append([data_by_organism['Run']
                                              .values[0]])
    #                 print(data_by_organism['MBases'])
            else:  # Check if there is transcriptomic data available
                if 'TRANSCRIPTOMIC' \
                        in data_by_organism.LibrarySource.values[0]:
                    sra_dic[organism] = [name_to_id[organism]]
                    data_o_t = data_by_organism[df.LibrarySource
                                                == 'TRANSCRIPTOMIC']
                    x = data_by_organism.iloc[(data_by_organism['MBases']-1000)
                                              .abs().argsort()[:1]]
                    max_cov_sample = x.MBases.values[0]
    #                 print(min_cov)
                    if max_cov_sample > 1000:  # for sure select the ones that have more than 10X coverage assuming a 1GB genome
                        all_index.append([x.index[0]])
                        sra_dic[organism].append(x['MBases'].values[0])
                        sra_dic[organism].append([x['Run'].values[0]])
                    else:
                        tmp = []
                        tmp_sra = []
                        cov = 0
                        for index, row in data_by_organism.iterrows():
                            cov += row['MBases']
                            if cov < 10001:
                                tmp.append(index)
                                tmp_sra.append(row['Run'])
                            else:
                                tmp_sra.append(row['Run'])
                                break
                        all_index.append(tmp)
                        sra_dic[organism].append(cov)
                        sra_dic[organism].append(tmp_sra)
                elif 'GENOMIC' in data_by_organism.LibrarySource.values[0]:
                    sra_dic[organism] = [name_to_id[organism]]
                    x = data_by_organism.iloc[(
                        data_by_organism['MBases']-10000).abs().argsort()[:1]]
                    x_all = data_by_organism.iloc[(
                        data_by_organism['MBases']-10000).abs().argsort()]
                    max_cov_sample = x.MBases.values[0]
                    if max_cov_sample > 10000:  # for sure select the ones that have more than 10X coverage assuming a 1GB genome
                        all_index.append([x.index[0]])
                        sra_dic[organism].append(x['MBases'].values[0])
                        sra_dic[organism].append([x['Run'].values[0]])
                    else:  # if there is not a single sample that fullfills coverage select mutliple runs
                        tmp = []
                        tmp_sra = []
                        cov = 0
                        for index, row in data_by_organism.iterrows():
                            cov += row['MBases']
                            if cov < 10001:
                                tmp.append(index)
                                tmp_sra.append(row['Run'])
                            else:
                                tmp_sra.append(row['Run'])
                                break
                        all_index.append(tmp)
                        sra_dic[organism].append(cov)
                        sra_dic[organism].append(tmp_sra)
    # return ordered dictionary by size of sra download
    sra_dic_ordered = OrderedDict(sorted(sra_dic.items(),
                                         key=lambda t: t[1][1]))
    if end == 0:  # use whole dictionary
        end = len(sra_dic_ordered.keys())
    return itertools.islice(sra_dic_ordered.items(), start, end)


def get_download_string(species_id, sra, se_pe='PAIRED'):
    sra_string = ''
    for i in sra:
        sra_string += '\"'+i+'\"'
        sra_string += ' '
    download = """#!/bin/bash
# $ -l mem=4G
# $ -S /bin/bash
# $ -l h_rt=10:00:0
# $ -pe smp 1
# $ -l tmpfs=100G
# $ -j y
# $ -N down_%s
# $ -wd /home/ucbpdvd/Scratch/avian/sge_output/
speciesid=%s
source activate r2t
mkdir /home/ucbpdvd/Scratch/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /home/ucbpdvd/Scratch/avian/reads/$speciesid
declare -a sra_all=(%s)
for sra in "${sra_all}"
do
    if [ "${sra:0:3}" == "SRR" ]; then
        ~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${sra:0:6}/$sra/$sra.sra ./
        echo 'Finished download'
        fastq-dump --split-files --gzip $sra.sra
        echo 'Finished fastq-dump'
        rm $sra.sra
    else
        ~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${sra:0:6}/002/$sra/$sra\_1.fastq.gz .
        ~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${sra:0:6}/002/$sra/$sra\_2.fastq.gz .
        echo 'Finished download'
    fi
done
cat *\_1.* > $speciesid\_1.fq.gz
cat *\_2.* > $speciesid\_2.fq.gz
rm SRR*
rm ERR*
echo 'Finished moving files'""" % (species_id, species_id, sra_string.rstrip())
    text_file = open('down_py_script.sh', "w")
    text_file.write(download)
    text_file.close()

    return 'down_py_script.sh'


def get_r2t_string(species_id, reference, se_pe='PAIRED', read_type='short'):
    if se_pe is 'PAIRED' and read_type is 'short':
        job_string = """#!/bin/bash
# $ -l mem=4G
# $ -S /bin/bash
# $ -l h_rt=12:00:0
# $ -pe smp 4
# $ -l tmpfs=120G
# $ -j y
# $ -N r2t_{species_id}
# $ -wd /home/ucbpdvd/Scratch/avian/sge_output/
reads=/home/ucbpdvd/Scratch/avian/reads/{species_id}
cd /home/ucbpdvd/Scratch/avian/r2t/
source activate r2t
python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fq.gz $reads/{species_id}_2.fq.gz \
--output_path /home/ucbpdvd/Scratch/avian/r2t/ --single_mapping {reference} \
--threads 4 --min_species 8 --read_type short --sample_reads \
--genome_len 1500000000 --coverage 10""".format(species_id=species_id,
                                                reference=reference)
    elif se_pe is 'SINGLE' and read_type is 'short':
        job_string = """#!/bin/bash
# $ -l mem=4G
# $ -S /bin/bash
# $ -l h_rt=12:00:0
# $ -pe smp 4
# $ -l tmpfs=120G
# $ -j y
# $ -N r2t_{species_id}
# $ -wd /home/ucbpdvd/Scratch/avian/sge_output/
reads=/home/ucbpdvd/Scratch/avian/reads/{species_id}
cd /home/ucbpdvd/Scratch/avian/r2t/
source activate r2t
python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree  \
--reads $reads/{species_id}_1.fq \
--output_path /home/ucbpdvd/Scratch/avian/r2t/ \
--single_mapping {reference} --threads 4 --min_species 8 --sample_reads \
--genome_len 1500000000 --coverage 10""".format(species_id=species_id,
                                                reference=reference)
    elif se_pe is 'SINGLE' and read_type is 'long':
        job_string = """#!/bin/bash
# $ -l mem=4G
# $ -S /bin/bash
# $ -l h_rt=12:00:0
# $ -pe smp 4
# $ -l tmpfs=120G
# $ -j y
# $ -N r2t_{species_id}
# $ -wd /home/ucbpdvd/Scratch/avian/sge_output/
reads=/home/ucbpdvd/Scratch/avian/reads/{species_id}
cd /home/ucbpdvd/Scratch/avian/r2t/
source activate r2t
python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fq \
--output_path /home/ucbpdvd/Scratch/avian/r2t/ \
--single_mapping {reference} --threads 4 --min_species 8 --read_type long \
--split_reads --sample_reads --genome_len 1500000000 \
--coverage 10""".format(species_id=species_id, reference=reference)

    text_file = open('r2t_py_script.sh', "w")
    text_file.write(job_string)
    text_file.close()

    return 'r2t_py_script.sh'


def get_rm_string(species_id):
    rm = """#!/bin/bash
# $ -l mem=4G
# $ -S /bin/bash
# $ -l h_rt=0:10:0
# $ -pe smp 1
# $ -j y
# $ -N rm_{species_id}
# $ -wd /home/ucbpdvd/Scratch/avian/sge_output/
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


def is_species_mapped(species_id, output):
    if os.path.exists(os.path.join(output, '03_mapping_'+species_id+'_1')):
        mapped_folder = os.path.join(output, '03_mapping_'+species_id+'_1')
        files = [f for f in glob.glob(os.path.join(mapped_folder, "*.fa"))]
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


def run_sge(sra_dic, output_speciesid):
    num_job_cycles = 0
    rm_job_id_idx = 0
    rm_job_id = []
    for species, sra in sra_dic:
        species_id = sra[0]
        sra_ids = sra[-1]
        # check whether the mapping already exists
        if not is_species_mapped(species_id, output_speciesid):
            # check whether the mapping already exists
            print('Submitting species {} with species id {}!'
                  .format(species, species_id))
            if num_job_cycles < 5:  # only run three jobs, then submit the jobs with dependency that files are again deleted
                # Set up download string
                job_string = get_download_string(species_id, sra_ids,
                                                 se_pe='PAIRED')

                # Open a pipe to the qsub command.
                p_download = output_shell('qsub ' + job_string)
                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = p_download.decode("utf-8").split(" ")[2]

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_speciesid,
                                                  '02_ref_dna/*.fa')):
                    # Set up r2t string
                    r2t_job_string = get_r2t_string(
                        species_id, ref, se_pe='PAIRED', read_type=sra[2])

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
                job_string = get_download_string(species_id, sra_ids,
                                                 se_pe='PAIRED')

                # Open a pipe to the qsub command.
                p_download = output_shell(
                    'qsub -hold_jid {} {}'.format(rm_job_id[rm_job_id_idx],
                                                  job_string))

                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = p_download.decode("utf-8").split(" ")[2]

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_speciesid,
                                                  '02_ref_dna/*.fa')):
                    # Set up r2t string
                    r2t_job_string = get_r2t_string(
                        species_id, ref, se_pe='PAIRED', read_type=sra[2])

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


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:o:a:e:",
                                   ["sra_file=", "out_speciesid=",
                                    "start_position=", "end_position="])
    except getopt.GetoptError as e:
        print(str(e))
        print('sge_submit.py -s <sra_file> -o <out_speciesid> '
              '-a <start_position> -e <end_position>')
        sys.exit(2)

    sra_file = None
    out_speciesid = None
    start = 0
    end = 0

    for opt, arg in opts:
        if opt == '-h':
            print('sge_submit.py -s <sra_file> -o <out_speciesid> '
                  '-a <start_position> -e <end_position>')
            sys.exit()
        elif opt in ("-s", "--sra_file"):
            sra_file = arg
        elif opt in ("-o", "--out_speciesid"):
            out_speciesid = arg
        elif opt in ("-a", "--start_position"):
            start = arg
        elif opt in ("-e", "--end_position"):
            end = arg
        else:
            assert False, "unhandled option"

    df = pd.read_csv(sra_file, sep='\t')
    df_illumina_paired = df[(df.Platform == 'ILLUMINA')
                            & (df.LibraryLayout == 'PAIRED')]
    print('Make sure to set the folders of reads, '
          'cluster ouput and location of run correctly!')

    print('We are selecting {} / {} illumina paired entries.'
          .format(len(set(df_illumina_paired.Organism)), len(df))

    # df = pd.read_csv(sra_file, sep='\t')
    name_to_id=get_name_to_id(df_illumina_paired)
    sra_dic=get_sra_dic(df_illumina_paired, name_to_id, start=start, end=end)

    run_sge(sra_dic, out_speciesid)


if __name__ == "__main__":
    main()
