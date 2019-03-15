from Bio import AlignIO
import pandas as pd
import ete3, os, shutil, glob, subprocess, tqdm
import logging.config, yaml
from pkg_resources import resource_string
import tempfile
from Bio.Align import MultipleSeqAlignment


from collections import Counter
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity

logging.getLogger(__name__).addHandler(logging.NullHandler())
conf = resource_string(__name__, '../read2tree/logging/log.yaml')

D = yaml.load(conf)
D.setdefault('version', 1)
logging.config.dictConfig(D)

logger = logging.getLogger(__name__)


def get_jaccard_sim(str1, str2):
    a = set(str1)
    b = set(str2)
    c = a.intersection(b)
    return float(len(c)) / (len(a) + len(b) - len(c))


def get_cosine_sim(*strs):
    vectors = [t for t in get_vectors(*strs)]
    return cosine_similarity(vectors)


def get_vectors(*strs):
    text = [t for t in strs]
    vectorizer = CountVectorizer(text)
    vectorizer.fit(text)
    return vectorizer.transform(text).toarray()

def get_alignment_sim(R2T_FOLDER, OUT_FOLDER, expected_nr_files=None):
    j_dist = {}
    cos_dist = {}
    with_r2t = {}
    num_sp_align_ref = {}
    num_sp_align_r2t = {}
    diff_num_sp = {}
    num_added = {}
    align_f1 = os.path.join(R2T_FOLDER, '03_align_aa')
    align_f2 = os.path.join(R2T_FOLDER, '06_align_merge_aa')
    if expected_nr_files:
        if len(glob.glob(os.path.join(align_f1, '*phy'))) == expected_nr_files:
            logger.info('All files for 03 are present!')
        else:
            f = open(os.path.join(OUT_FOLDER, os.path.abspath(R2T_FOLDER).split("/")[-1] + 'ERROR.csv'), "w+")
            f.write('that did not work')
        if len(glob.glob(os.path.join(align_f2, '*fa'))) == expected_nr_files:
            logger.info('All files for 06 are present!')
        else:
            f = open(os.path.join(OUT_FOLDER, os.path.abspath(R2T_FOLDER).split("/")[-1] + 'ERROR.csv'), "w+")
            f.write('that did not work')

    for a1 in glob.glob(os.path.join(align_f1, '*phy')):
        og = os.path.basename(a1).split('.')[0]
        align1 = AlignIO.read(a1, 'phylip-relaxed')
        a2 = os.path.join(align_f2, og + '.fa')
        align2 = AlignIO.read(a2, 'phylip-relaxed')
        with_r2t[og] = "|".join([r.id for r in align2 if '_1' in r.id])
        num_added[og] = len([r.id for r in align2 if '_1' in r.id])
        num_sp_align_ref[og] = len(align1)
        num_sp_align_r2t[og] = len(align2)
        diff_num_sp[og] = len(align1) - len(align2)
        if len(align1) == len(align2):
            align1.sort()
            align2.sort()
            align1_words = [align1[:, i] for i in range(len(align1[1, :]))]
            align2_words = [align2[:, i] for i in range(len(align2[1, :]))]
            j_dist[og] = get_jaccard_sim(align1_words, align2_words)
            cos_dist[og] = get_cosine_sim(' '.join(align1_words), ' '.join(align2_words))[-1,0]
        elif len(align2) < len(align1):  # in case the alingments are not filled
            sp_align2 = [r.id[0:5] for r in align2]
            new_align1 = MultipleSeqAlignment([r for r in align1 if r.id in sp_align2])
            new_align1.sort()
            align2.sort()
            new_align1_words = [new_align1[:, i] for i in range(len(new_align1[1, :]))]
            align2_words = [align2[:, i] for i in range(len(align2[1, :]))]
            j_dist[og] = get_jaccard_sim(new_align1_words, align2_words)
            cos_dist[og] = get_cosine_sim(' '.join(new_align1_words), ' '.join(align2_words))[-1, 0]
        elif len(align2) > len(align1):  # in case the alingments are not filled
            sp_align1 = [r.id[0:5] for r in align1]
            new_align2 = MultipleSeqAlignment([r for r in align2 if r.id[0:5] in sp_align1])
            new_align2.sort()
            align1.sort()
            new_align2_words = [new_align2[:, i] for i in range(len(new_align2[1, :]))]
            align1_words = [align1[:, i] for i in range(len(align1[1, :]))]
            j_dist[og] = get_jaccard_sim(new_align2_words, align1_words)
            cos_dist[og] = get_cosine_sim(' '.join(align1_words), ' '.join(new_align2_words))[-1, 0]
    df = pd.DataFrame({'num_added': pd.Series(num_added), 'cos_dist':pd.Series(cos_dist), 'j_dist': pd.Series(j_dist), 'with_r2t': pd.Series(with_r2t),
                       'num_sp_align_ref': pd.Series(num_sp_align_ref), 'num_sp_align_r2t': pd.Series(num_sp_align_r2t),
                       'diff_num_sp': pd.Series(diff_num_sp)})
    df.to_csv(
        os.path.join(OUT_FOLDER, os.path.abspath(R2T_FOLDER).split("/")[-1] + '_alignments.csv'))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""Compute column alignment distances """)

    parser.add_argument('--r2t_folder', default=None,
                        help='r2t_folder')
    parser.add_argument('--output_folder', default=None,
                        help='output_file')
    parser.add_argument('--expected_nr_files', default=None,
                        help='expected_nr_files')

    conf = parser.parse_args()
    logger.info('----- Lets get the alignment scores -----')
    get_alignment_sim(conf.r2t_folder, conf.output_folder, expected_nr_files=int(conf.expected_nr_files))
