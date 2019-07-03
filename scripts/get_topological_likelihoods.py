from Bio import AlignIO
import pandas as pd
import ete3, os, shutil, glob, subprocess, tqdm
import logging.config, yaml
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
        logger.info("Shell command failed to execute")
        logger.info(line)
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

def collect_likelihood(tmp_folder, suffix, method='raxml'):
    if 'raxml' in method:
        run = os.path.join(tmp_folder, 'RAxML_info.'+suffix)
        with open(run, 'r') as f:
            tmp = f.readlines()
        return float([l.split(' ')[-1].replace('\n','') for l in tmp if 'Final' in l][0])
    else:
        run = os.path.join(tmp_folder, suffix+'.log')
        with open(run, 'r') as f:
            tmp = f.readlines()
        try:
            value = float([l.split(' ')[-1].replace('\n','') for l in tmp if 'BEST SCORE FOUND' in l][0])
        except IndexError:
            print(tmp)
            value = 0.0
        return value

def get_likelihoods(R2T_FOLDER, TREE_FILE1, TREE_FILE2, OUT_FOLDER, tree_method='iqtree'):

    likelihood_refa_reft = {}
    likelihood_refa_r2tt = {}
    likelihood_diff_refa = {}
    likelihood_r2ta_reft = {}
    likelihood_r2ta_r2ta = {}
    likelihood_diff_r2ta = {}

    align_f1 = os.path.join(R2T_FOLDER, '03_align_aa')
    align_f2 = os.path.join(R2T_FOLDER, '06_align_merge_aa')

    for f in tqdm.tqdm(glob.glob(os.path.join(align_f1, '*.phy'))):
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
        #         logger.info(os.path.getsize(fixed_align_out))

        trim_tree(TREE_FILE1, species, tmp_tree1_out)
        if 'raxml' in tree_method:
            output_shell(
                'raxmlHPC -g ' + tmp_tree1_out + ' -p 12345 -T 4 -m PROTGAMMALGX -s ' + fixed_align_out + ' -n rax1 -w' + tmp_folder)
            likelihood_refa_reft[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, fixed_align_out, method='raxml')
        elif 'iqtree' in tree_method:
            output_shell(
                'iqtree -te ' + tmp_tree1_out + ' -redo -seed 12345 -nt 4 -mem 4G -seed 12345 -m LG -s ' + fixed_align_out)
            likelihood_refa_reft[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, fixed_align_out,
                                                                                         method='iqtree')

        trim_tree(TREE_FILE2, species, tmp_tree2_out)
        if 'raxml' in tree_method:
            output_shell(
                'raxmlHPC -g ' + tmp_tree2_out + ' -p 12345 -T 4 -m PROTGAMMALGX -s ' + fixed_align_out + ' -n rax2 -w' + tmp_folder)
            likelihood_refa_r2tt[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, fixed_align_out, method='raxml')
        elif 'iqtree' in tree_method:
            output_shell(
                'iqtree -te ' + tmp_tree2_out + ' -redo -seed 12345 -nt 4 -mem 4G -seed 12345 -m LG -s ' + fixed_align_out)
            likelihood_refa_r2tt[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, fixed_align_out,
                                                                                         method='iqtree')
        likelihood_diff_refa[os.path.basename(f).split('.')[0]] = likelihood_refa_reft[
                                                                          os.path.basename(f).split('.')[0]] - \
                                                                      likelihood_refa_r2tt[
                                                                          os.path.basename(f).split('.')[0]]

        shutil.rmtree(tmp_folder)

    for f in tqdm.tqdm(glob.glob(os.path.join(align_f2, '*.fa'))):
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
        if 'raxml' in tree_method:
            output_shell(
                'raxmlHPC -g ' + tmp_tree1_out + ' -p 12345 -T 4 -m PROTGAMMALGX -s ' + fixed_align_out + ' -n rax1 -w' + tmp_folder)
            likelihood_r2ta_reft[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, fixed_align_out, method='raxml')
        elif 'iqtree' in tree_method:
            output_shell(
                'iqtree -te ' + tmp_tree1_out + ' -redo -seed 12345 -nt 4 -mem 4G -seed 12345 -m LG -s ' + fixed_align_out)
            likelihood_r2ta_reft[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, fixed_align_out,
                                                                                         method='iqtree')

        trim_tree(TREE_FILE2, species, tmp_tree2_out)
        if 'raxml' in tree_method:
            output_shell(
                'raxmlHPC -g ' + tmp_tree2_out + ' -p 12345 -T 4 -m PROTGAMMALGX -s ' + fixed_align_out + ' -n rax2 -w' + tmp_folder)
            likelihood_r2ta_r2ta[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, fixed_align_out, method='raxml')
        elif 'iqtree' in tree_method:
            output_shell(
                'iqtree -te ' + tmp_tree2_out + ' -redo -seed 12345 -nt 4 -mem 4G -seed 12345 -m LG -s ' + fixed_align_out)
            likelihood_r2ta_r2ta[os.path.basename(f).split('.')[0]] = collect_likelihood(tmp_folder, fixed_align_out,
                                                                                         method='iqtree')

        likelihood_diff_r2ta[os.path.basename(f).split('.')[0]] = likelihood_r2ta_reft[
                                                                          os.path.basename(f).split('.')[0]] - \
                                                                      likelihood_r2ta_r2ta[
                                                                          os.path.basename(f).split('.')[0]]

        shutil.rmtree(tmp_folder)

    df = pd.DataFrame({'likelihood_refa_reft': pd.Series(likelihood_refa_reft),
                       'likelihood_refa_r2tt': pd.Series(likelihood_refa_r2tt),
                       'likelihood_diff_refa': pd.Series(likelihood_diff_refa),
                       'likelihood_r2ta_reft': pd.Series(likelihood_r2ta_reft),
                       'likelihood_r2ta_r2ta': pd.Series(likelihood_r2ta_r2ta),
                       'likelihood_diff_r2ta': pd.Series(likelihood_diff_r2ta)})
    df.to_csv(
        os.path.join(OUT_FOLDER, os.path.abspath(R2T_FOLDER).split("/")[-1] + '_likelihoods.csv'))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""Compute raxml likelihood scores for """)

    parser.add_argument('--ref_tree', default=None,
                        help='ref_tree')
    parser.add_argument('--r2t_tree', default=None,
                        help='r2t_tree')
    parser.add_argument('--r2t_folder', default=None,
                        help='r2t_folder')
    parser.add_argument('--tree_method', default=None,
                        help='tree_method')
    parser.add_argument('--output_folder', default=None,
                        help='output_file')

    conf = parser.parse_args()
    logger.info('----- Lets get the likelihoods -----')
    get_likelihoods(conf.r2t_folder, conf.ref_tree, conf.r2t_tree, conf.output_folder, tree_method=conf.tree_method)
