from Bio import AlignIO
import pandas as pd
import ete3, os, shutil, glob, subprocess, tqdm, tempfile
import logging, yaml
from pkg_resources import resource_string
import tempfile

logging.getLogger(__name__).addHandler(logging.NullHandler())
conf = resource_string(__name__, '../read2tree/logging/log.yaml')

D = yaml.load(conf)
D.setdefault('version', 1)
logging.config.dictConfig(D)

logger = logging.getLogger(__name__)

def output_shell(line):
    """
    Save output of shell line that has pipes
    taken from: https://stackoverflow.com/questions/7389662/link-several-popen-commands-with-pipes
    :param line:
    :return:
    """
    try:
        shell_command = subprocess.Popen(
            line, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            shell=True)
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

def trim_tree(t_file, species, tree_out):
    f = open(t_file, 'r')
    t1 = ete3.Tree(f.read())
    for leaf in t1:
        if '_1' in leaf.name:
            tmp = leaf.name[0:5]
            leaf.name = tmp
    t1.prune(species)
    t1.write(format=1, outfile=tree_out)

def collect_likelihood(tmp_folder, suffix):
    raxml_run = os.path.join(tmp_folder, 'RAxML_info.'+suffix)
    with open(raxml_run, 'r') as f:
        tmp = f.readlines()
    return float([l.split(' ')[-1].replace('\n','') for l in tmp if 'Final' in l][0])

def get_likelihoods(ALIGN_FOLDER, TREE_FILE1, TREE_FILE2, outfile):
    t1_dic = {}
    t2_dic = {}
    for f in tqdm.tqdm(glob.glob(os.path.join(ALIGN_FOLDER, '*.fa'))):
        tmp_folder = tempfile.mkdtemp()

        fixed_align_out = os.path.join(tmp_folder, os.path.basename(f))
        tmp_tree1_out = os.path.join(tmp_folder, 't1.nwk')
        tmp_tree2_out = os.path.join(tmp_folder, 't2.nwk')

        align = AlignIO.read(f, 'phylip-relaxed')
        species = []
        for r in align:
            tmp = r.id[0:5]
            r.id = tmp
            species.append(tmp)

        AlignIO.write(align, fixed_align_out, 'phylip')
        #     if os.path.getsize(fixed_align_out) > 0:
        #         print(os.path.getsize(fixed_align_out))

        trim_tree(TREE_FILE1, species, tmp_tree1_out)
        output_shell(
            'raxmlHPC -g ' + tmp_tree1_out + ' -p 12345 -m PROTGAMMALGX -s ' + fixed_align_out + ' -n rax1 -w' + tmp_folder)
        t1_dic[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, 'rax1')

        trim_tree(TREE_FILE2, species, tmp_tree2_out)
        output_shell(
            'raxmlHPC -g ' + tmp_tree2_out + ' -p 12345 -m PROTGAMMALGX -s ' + fixed_align_out + ' -n rax2 -w' + tmp_folder)
        t2_dic[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, 'rax2')

        shutil.rmtree(tmp_folder)

    ogs = list(t1_dic.keys())
    t1_likelihoods = [t1_dic[k] for k in og_s]
    t2_likelihoods = [t2_dic[k] for k in og_s]

    df = pd.Dataframe({'OGs': ogs, 'T1_likelihoods':t1_likelihoods, 'T2_likelihoods':t2_likelihoods})
    df.to_csv(outfile)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""Compute raxml likelihood scores for """)

    parser.add_argument('--tree1', default=None,
                        help='tree1')
    parser.add_argument('--tree2', default=None,
                        help='tree2')
    parser.add_argument('--alignment_folder', default=None,
                        help='alignment_folder)
    parser.add_argument('--output_file', default=None,
                        help='output_file')

    conf = parser.parse_args()
    logger.info('----- WE ARE SAMPLING THE READS -----')
    get_likelihoods(conf.alignment_folder, conf.tree1, conf.tree2, conf.output_file)
