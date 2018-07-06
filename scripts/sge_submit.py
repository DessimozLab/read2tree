import os
import sys
import getopt
import glob
import time
import subprocess
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
                    x = data_by_organism.iloc[(
                        data_by_organims['MBases'] - 1000).abs().argsort()[:1]]
                    if x.MBases.values[
                            0] > 1000:  # for sure select the ones that have more than 10X coverage assuming a 1GB genome
                        all_index.append(x.index[0])
                elif 'GENOMIC' in data_by_organism.LibrarySource:
                    x = data_by_organism.iloc[(
                        data_by_organims['MBases'] - 10000).abs().argsort()[:1]]
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
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=10:00:0
#$ -pe smp 1
#$ -l tmpfs=100G
#$ -j y
#$ -N down_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
srr={sra}
speciesid={species_id}
source activate r2t
mkdir /home/ucbpdvd/Scratch/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{sra_p}/002/$srr/$srr\_1.fastq.gz .
~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{sra_p}/002/$srr/$srr\_2.fastq.gz .
echo 'Finished download'
mv $srr\_1.fastq.gz $speciesid\_1.fq.gz
mv $srr\_2.fastq.gz $speciesid\_2.fq.gz
echo 'Finished moving files'""".format(species_id=species_id, sra=sra,
                                       sra_p=sra[0:6])
    if 'ERR' in sra and se_pe is 'SINGLE':
        download = """#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=10:00:0
#$ -pe smp 1
#$ -l tmpfs=100G
#$ -j y
#$ -N down_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
srr={sra}
speciesid={species_id}
source activate r2t
mkdir /home/ucbpdvd/Scratch/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /scratch/beegfs/weekly/ddylus/avian/reads/$speciesid
~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{sra_p}/002/$srr/$srr.fastq.gz .
echo 'Finished download'
mv $srr.fastq.gz $speciesid\_1.fq.gz
echo 'Finished moving files'""".format(species_id=species_id, sra=sra,
                                       sra_p=sra[0:6])
    elif 'SRR' in sra and se_pe is 'PAIRED':
        download = """#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=10:00:0
#$ -pe smp 4
#$ -l tmpfs=100G
#$ -j y
#$ -N down_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
srr={sra}
speciesid={species_id}
source activate r2t
mkdir /home/ucbpdvd/Scratch/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /home/ucbpdvd/Scratch/avian/reads/$speciesid
~/.aspera/connect/bin/ascp  -v -QT -k1 -l100M -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/{sra_p}/$srr/$srr.sra ./
echo 'Finished download'
fastq-dump --split-files --gzip $srr.sra
echo 'Finished getting fastq from sra and split files'
rm *.sra
mv *\_1.* $speciesid\_1.fq.gz
mv *\_2.* $speciesid\_2.fq.gz
echo 'Finished moving files'""".format(species_id=species_id, sra=sra,
                                       sra_p=sra[0:6])

    elif 'SRR' in sra and se_pe is 'SINGLE':
        download = """#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=10:00:0
#$ -pe smp 4
#$ -l tmpfs=100G
#$ -j y
#$ -N down_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
srr={sra}
speciesid={species_id}
source activate r2t
mkdir /home/ucbpdvd/Scratch/avian/reads/$speciesid
echo 'Created read $speciesid'
cd /home/ucbpdvd/Scratch/avian/reads/$speciesid
~/.aspera/connect/bin/ascp -v -QT -k1 -l100M -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/{sra_p}/$srr/$srr.sra ./
echo 'Finished download'
fastq-dump --gzip *.sra
rm *.sra
echo 'Finished getting fastq from sra'
mv *.gz $speciesid\_1.fq.gz
echo 'Finished moving files'""".format(species_id=species_id, sra=sra,
                                       sra_p=sra[0:6])

    text_file = open('down_py_script.sh', "w")
    text_file.write(download)
    text_file.close()

    return 'down_py_script.sh'


def get_r2t_string(species_id, reference, se_pe='PAIRED', read_type='short'):
    if se_pe is 'PAIRED' and read_type is 'short':
        job_string = '''#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=12:00:0
#$ -pe smp 4
#$ -l tmpfs=120G
#$ -j y
#$ -N r2t_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
reads=/home/ucbpdvd/Scratch/avian/reads/{species_id}
cd /home/ucbpdvd/Scratch/avian/r2t/
source activate r2t
python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fq.gz $reads/{species_id}_2.fq.gz \
--output_path /home/ucbpdvd/Scratch/avian/r2t/ --single_mapping {reference} \
--threads 4 --min_species 8 --read_type short --sample_reads \
--genome_len 1500000000 --coverage 10'''.format(species_id=species_id,
                                                reference=reference)
    elif se_pe is 'SINGLE' and read_type is 'short':
        job_string = '''#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=12:00:0
#$ -pe smp 4
#$ -l tmpfs=120G
#$ -j y
#$ -N r2t_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
reads=/home/ucbpdvd/Scratch/avian/reads/{species_id}
cd /home/ucbpdvd/Scratch/avian/r2t/
source activate r2t
python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree  \
--reads $reads/{species_id}_1.fq \
--output_path /home/ucbpdvd/Scratch/avian/r2t/ \
--single_mapping {reference} --threads 4 --min_species 8 --sample_reads \
--genome_len 1500000000 --coverage 10'''.format(species_id=species_id,
                                                reference=reference)
    elif se_pe is 'SINGLE' and read_type is 'long':
        job_string = '''#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=12:00:0
#$ -pe smp 4
#$ -l tmpfs=120G
#$ -j y
#$ -N r2t_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
reads=/home/ucbpdvd/Scratch/avian/reads/{species_id}
cd /home/ucbpdvd/Scratch/avian/r2t/
source activate r2t
python -W ignore /home/ucbpdvd/opt/read2tree/bin/read2tree \
--reads $reads/{species_id}_1.fq \
--output_path /home/ucbpdvd/Scratch/avian/r2t/ \
--single_mapping {reference} --threads 4 --min_species 8 --read_type long \
--split_reads --sample_reads --genome_len 1500000000 \
--coverage 10'''.format(species_id=species_id, reference=reference)

    text_file = open('r2t_py_script.sh', "w")
    text_file.write(job_string)
    text_file.close()

    return 'r2t_py_script.sh'


def get_rm_string(species_id):
    rm = '''#!/bin/bash
#$ -l mem=4G
#$ -S /bin/bash
#$ -l h_rt=0:10:0
#$ -pe smp 1
#$ -j y
#$ -N rm_{species_id}
#$ -wd /home/ucbpdvd/Scratch/avian/sge_output/
rm -r \
/home/ucbpdvd/Scratch/avian/reads/{species_id}'''.format(species_id=species_id)

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


def run_sge(sra_dic, output_speciesid):
    num_job_cycles = 0
    rm_job_id_idx = 0
    rm_job_id = []
    for species, sra in sra_dic.items():
        species_id = sra[-1]
        # check whether the mapping already exists
        if not is_species_mapped(species_id, output_speciesid):
            # check whether the mapping already exists
            print('Submitting species {} with species id {}!'.format(species, species_id))
            if num_job_cycles < 3:  # only run three jobs, then submit the jobs with dependency that files are again deleted
                # Set up download string
                job_string = get_download_string(species_id, sra[0], se_pe=sra[1])

                # Open a pipe to the qsub command.
                p_download = output_shell('qsub ' + job_string)
                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = p_download.decode("utf-8").split(" ")[2]

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_speciesid, '02_ref_dna/*.fa')):
                    # Set up r2t string
                    r2t_job_string = get_r2t_string(species_id, ref, se_pe=sra[1], read_type=sra[2])

                    # Open a pipe to the qsub command.
                    #output_r2t, input_r2t = Popen('qsub -hold_jid {}'.format(jobid))
                    p_download = output_shell('qsub -hold_jid {} {}'.format(jobid, r2t_job_string))

                    # Append jobid of r2t
                    r2t_jobids.append(p_download.decode("utf-8").split(" ")[2])

                    time.sleep(0.1)

                # Set up r2t string
                rm_job_string = get_rm_string(species_id)

                # Open a pipe to the qsub command.
                #output_rm, input_rm = Popen('qsub -hold_jid {}'.format(','.join(r2t_jobids)))
                p_rm = output_shell(
                    'qsub -hold_jid {} {}'.format(','.join(r2t_jobids), rm_job_string))

                # Print your job and the system response to the screen as it's submitted
                rm_job_id.append(p_rm.decode("utf-8").split(" ")[2])
                time.sleep(0.1)
            else:  # this part should ensure that on the scratch never more than 3 downloads are available
                # Set up download string
                job_string = get_download_string(species_id, sra[0], se_pe=sra[1])

                # Open a pipe to the qsub command.
                p_download = output_shell(
                    'qsub -hold_jid {} {}'.format(rm_job_id[rm_job_id_idx], job_string))

                time.sleep(0.1)

                # Get jobid for job chaining
                jobid = p_download.decode("utf-8").split(" ")[2]

                r2t_jobids = []
                for ref in glob.glob(os.path.join(output_speciesid, '02_ref_dna/*.fa')):
                    # Set up r2t string
                    r2t_job_string = get_r2t_string(species_id, ref, se_pe=sra[1], read_type=sra[2])

                    # Open a pipe to the qsub command.
                    # output_r2t, input_r2t = Popen('qsub -hold_jid {}'.format(jobid))
                    p_download = output_shell('qsub -hold_jid {} {}'.format(jobid, r2t_job_string))

                    # Append jobid of r2t
                    r2t_jobids.append(p_download.decode("utf-8").split(" ")[2])

                    time.sleep(0.1)

                # Set up r2t string
                rm_job_string = get_rm_string(species_id)

                # Open a pipe to the qsub command.
                # output_rm, input_rm = Popen('qsub -hold_jid {}'.format(','.join(r2t_jobids)))
                p_rm = output_shell(
                    'qsub -hold_jid {} {}'.format(','.join(r2t_jobids), rm_job_string))

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
            print('sge_submit.py -s <sra_file> -o <out_speciesid>')
            sys.exit()
        elif opt in ("-s", "--sra_file"):
            sra_file = arg
        elif opt in ("-o", "--out_speciesid"):
            out_speciesid = arg
        else:
            assert False, "unhandled option"

    # df = pd.read_csv(sra_file, sep='\t')
    sra_dic = {'Zonotrichia leucophrys': ['SRR1199463', 'PAIRED', 'short', 'ZONLE'],
               'Periparus ater ater': ['SRR1810767', 'PAIRED', 'short', 'PERA0'],
               'Phylloscopus trochiloides trochiloides': ['SRR3217927',
                                                          'PAIRED',
                                                          'short',
                                                          'PHYTR'],
               'Periparus ater sardus': ['SRR1810764', 'PAIRED', 'short', 'PERA1'],
               'Callipepla californica': ['SRR3630008', 'PAIRED', 'short', 'CALCA'],
               'Sibirionetta formosa': ['SRR3471611', 'PAIRED', 'short', 'SIBFO'],
               'Nesoptilotis leucotis': ['SRR3901721', 'PAIRED', 'short', 'NESLE'],
               'Macronectes giganteus': ['SRR6902605', 'PAIRED', 'short', 'MACGI'],
               'Cyanistes caeruleus calamensis': ['SRR1810758', 'PAIRED', 'short', 'CYACA'],
               'Sericornis frontalis': ['SRR3901710', 'PAIRED', 'short', 'SERFR'],
               'Cyanistes palmensis': ['SRR1810654', 'PAIRED', 'short', 'CYAPA'],
               'Periparus ater atlas': ['SRR1810770', 'PAIRED', 'short', 'PERA2'],
               'Falco cherrug cherrug': ['SRR671935', 'SINGLE', 'short', 'FALCH'],
               'Akialoa obscura': ['SRR3181052', 'SINGLE', 'short', 'AKIOB'],
               'Phylloscopus trochiloides viridanus': ['SRR1172475',
                                                       'PAIRED',
                                                       'short',
                                                       'PHYT0'],
               'Ammodramus caudacutus': ['SRR1957204', 'PAIRED', 'short', 'AMMCA'],
               'Myadestes myadestinus': ['SRR3181004', 'SINGLE', 'short', 'MYAMY'],
               'Cyanistes flavipectus': ['SRR1810646', 'PAIRED', 'short', 'CYAFL'],
               'Callaeas cinereus': ['SRR3180908', 'SINGLE', 'short', 'CALCI'],
               'Phylloscopus plumbeitarsus': ['SRR3223375', 'PAIRED', 'short', 'PHYPL'],
               'Platycercus eximius': ['SRR3901724', 'PAIRED', 'short', 'PLAEX'],
               'Lonchura leucosticta': ['SRR5976562', 'PAIRED', 'short', 'LONLE'],
               'Bubo bubo': ['SRR3203225', 'PAIRED', 'short', 'BUBBU'],
               'Campylorhynchus brunneicapillus': ['SRR1145744', 'SINGLE', 'short', 'CAMBR'],
               'Ammodramus nelsoni': ['SRR1955654', 'PAIRED', 'short', 'AMMNE'],
               'Phasianidae gen. sp.': ['SRR088928', 'SINGLE', 'long', 'PHAGE'],
               'Pelecanus occidentalis': ['SRR1145758', 'PAIRED', 'short', 'PELOC'],
               'Zonotrichia leucophrys gambelii': ['SRR1238747', 'PAIRED', 'short', 'ZONL0'],
               'Athene noctua': ['SRR3203242', 'PAIRED', 'short', 'ATHNO'],
               'Gymnopithys rufigula': ['SRR3115006', 'PAIRED', 'short', 'GYMRU'],
               'Machlolophus spilonotus': ['SRR765718', 'PAIRED', 'short', 'MACSP'],
               'Zosterops lateralis': ['SRR2145253', 'SINGLE', 'short', 'ZOSLA'],
               'Cnemophilus loriae': ['SRR2968794', 'PAIRED', 'short', 'CNELO'],
               'Alca torda': ['SRR1145756', 'PAIRED', 'short', 'ALCTO'],
               'Accipiter virgatus': ['SRR3203234', 'PAIRED', 'short', 'ACCVI'],
               'Lanius excubitor': ['SRR2968729', 'PAIRED', 'short', 'LANEX'],
               'Paradoxornis webbianus bulomachus': ['SRR392516',
                                                     'SINGLE',
                                                     'short',
                                                     'PARWE'],
               'Cyanistes caeruleus caeruleus': ['SRR1810777', 'PAIRED', 'short', 'CYAC0'],
               'Butastur indicus': ['SRR3203233', 'PAIRED', 'short', 'BUTIN'],
               'Turnagra capensis': ['SRR3180914', 'SINGLE', 'short', 'TURCA'],
               'Mareca falcata': ['SRR3471610', 'PAIRED', 'short', 'MARFA'],
               'Nucifraga columbiana': ['SRR1166560', 'PAIRED', 'short', 'NUCCO'],
               'Otus scops': ['SRR3203230', 'PAIRED', 'short', 'OTUSC'],
               'Malurus lamberti': ['SRR3901709', 'PAIRED', 'short', 'MALLA'],
               'Passer montanus': ['SRR5369936', 'PAIRED', 'short', 'PASMO'],
               'Zimmerius gracilipes': ['SRR1021716', 'PAIRED', 'short', 'ZIMGR'],
               'Heteralocha acutirostris': ['SRR3180975', 'SINGLE', 'short', 'HETAC'],
               'Falco tinnunculus': ['SRR3203231', 'PAIRED', 'short', 'FALTI'],
               'Cyanistes caeruleus ogliastrae': ['SRR1810757', 'PAIRED', 'short', 'CYAC1'],
               'Picus canus': ['SRR3203240', 'PAIRED', 'short', 'PICCA'],
               'Falco peregrinus peregrinus': ['SRR671934', 'SINGLE', 'short', 'FALPE'],
               'Asio otus': ['SRR3203220', 'PAIRED', 'short', 'ASIOT'],
               'Anser sp.': ['SRR1060398', 'PAIRED', 'short', 'ANSSP'],
               'Spinus cucullatus': ['SRR2895762', 'PAIRED', 'short', 'SPICU'],
               'Psephotellus pulcherrimus': ['SRR3180905', 'SINGLE', 'short', 'PSEPU'],
               'Lagopus lagopus': ['SRR2913174', 'PAIRED', 'short', 'LAGLA'],
               'Glyphorynchus spirurus': ['SRR3115005', 'PAIRED', 'short', 'GLYSP']}

    run_sge(sra_dic, out_speciesid)


if __name__ == "__main__":
    main()
