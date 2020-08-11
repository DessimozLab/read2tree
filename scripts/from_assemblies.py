#!/usr/bin/env python3
import os

import Bio.AlignIO
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy
import collections
import logging
logger = logging.getLogger(__name__)


def chunk_msa(msa, min_length=300):
    def frac_uninform(column):
        cnts = collections.Counter(column)
        return (cnts['-'] + cnts['n']) / len(column)

    uninform = numpy.array([frac_uninform(msa[:, i]) for i in range(msa.get_alignment_length())])
    splits = []
    cur_chunk = 0
    for i in range(msa.get_alignment_length()):
        if uninform[i] <= 0.1:
            cur_chunk += 1
        if cur_chunk >= min_length:
            splits.append(i)
            cur_chunk = 0
    # remove last split such that the last segement is still at least min_length long
    return splits[:-1]


def split_msa(msa, split_pos):
    msas = []
    for i in range(len(split_pos)+1):
        rng = slice(split_pos[i-1] if i > 0 else None,
                    split_pos[i] if i < len(split_pos) else None)
        cur_msa_chunk = msa[:, rng]
        for rec in cur_msa_chunk:
            rec.id = rec.id + "_OG{}".format(i)
        msas.append(cur_msa_chunk)
    logger.info("split msa into {} chunks".format(len(msas)))
    assert sum((x.get_alignment_length() for x in msas)) == msa.get_alignment_length()
    return msas


def write_aligned_og(msas, outdir):
    os.makedirs(outdir)
    for og, msa in enumerate(msas):
        fname = os.path.join(outdir, "OG{}.phy".format(og))
        with open(fname, 'wt') as fout:
            for rec in msa:
                rec.id = rec.id.split('_OG')[0]
            Bio.AlignIO.write(msa, fout, 'phylip')


def get_unaligned_seq_rec_from_aligned_seq(rec):
    seq = Seq("".join([bp for bp in str(rec.seq) if bp != '-']),
              alphabet=rec.seq.alphabet)
    new = SeqRecord(seq)
    new.id = rec.id
    new.name = rec.name
    return new


def compute_seq_identity(s1, s2):
    assert(len(s1) == len(s2))
    x = 0
    for i, j in zip(s1, s2):
        if i == '-' and j == 'j':
            continue
        if i == j:
           x += 1
    return x / len(s1)


def select_representative_references(msa, min_divergence):
    if min_divergence <= 0:
        return [s.id for s in msa]
    rep = [0]
    for k, rec in enumerate(msa):
        if k == 0: continue
        for ref in rep:
            if 100 * compute_seq_identity(rec, msa[ref]) > 100-min_divergence:
                break
        else:
            rep.append(k)
    return [msa[i].id for i in rep]


def write_unaligned_og(msas, outdir):
    os.makedirs(outdir)
    for og, msa in enumerate(msas):
        fname = os.path.join(outdir, "OG{}.fa".format(og))
        with open(fname, 'wt') as fout:
            for rec in msa:
                Bio.SeqIO.write(get_unaligned_seq_rec_from_aligned_seq(rec), fout, 'fasta')


def write_reference_files(msas, outdir):
    os.makedirs(outdir)
    for rg in range(len(msas[0])):
        fname = os.path.join(outdir, "{}_OGs.fa".format(msas[0][rg].id.split('_')[0]))
        with open(fname, 'wt') as fout:
            for msa in msas:
                Bio.SeqIO.write(get_unaligned_seq_rec_from_aligned_seq(msa[rg]), fout, 'fasta')


def make_genome_names_compatible(msa, out):
    with open(os.path.join(out, "genome_mapping.txt"), 'wt') as fh:
        for i in range(len(msa)):
            orig_header = msa[i].description
            msa[i].description = ""
            msa[i].id = "RG{:03d}".format(i+1)
            fh.write("RG{:03d}\t{}\n".format(i+1, orig_header))
    return msa


def prepare_r2t(msa, base_path, min_divergence):
    os.makedirs(base_path)
    msa = make_genome_names_compatible(msa, base_path)
    rep_refs = select_representative_references(msa, min_divergence)
    split_pos = chunk_msa(msa)
    all_msas = split_msa(msa, split_pos)
    write_unaligned_og(all_msas, os.path.join(base_path, "01_ref_ogs_dna"))
    write_reference_files(all_msas, os.path.join(base_path, "02_ref_dna"))
    write_aligned_og(all_msas, os.path.join(base_path, "03_align_dna"))
    try:
        fd = os.open(base_path, os.O_RDONLY)
        os.symlink("01_ref_ogs_dna", "01_ref_ogs_aa", dir_fd=fd)
        os.symlink("03_align_dna", "03_align_aa", dir_fd=fd)
    finally:
        os.close(fd)

    with open(os.path.join(base_path, "mplog.log"), 'w') as log:
        log.write("2020-04-18 00:30:16,393 - read2tree.ReferenceSet - INFO - test_1a: "
                  "Extracted {} reference species form 6428 ogs took 0.045595407485961914\n".format(len(msa)))
        log.write("2020-04-18 00:30:16,338 - read2tree.OGSet - INFO - test_1a: "
                  "Gathering of DNA seq for {} OGs took 1037.3560745716095.\n".format(len(all_msas)))
        log.write("2020-04-18 01:38:40,309 - read2tree.Aligner - INFO - test_1a: "
                  "Alignment of {} OGs took 4103.362480401993.\n".format(len(all_msas)))
    print("We suggest to use the following reference genomes for mapping:")
    print(rep_refs)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="extract from a whole genome alignment MSA pseudo OGs for r2t")
    parser.add_argument('--out', required=True, help="folder where output is written")
    parser.add_argument('-d', '--min-divergence', type=int, default=0,
                        help="minimum divergence (seq-identity) for reference genome to be used for "
                             "separate mapping. [Default 0]")
    parser.add_argument('-v', action='count', help="increase verbosity")
    parser.add_argument('msa', help="input whole genome alignment MSA in fasta/phylip format")
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))

    if conf.msa.endswith('.fa') or conf.msa.endswith('.fasta'):
        format = 'fasta'
    elif conf.msa.endswith('.phy'):
        format = 'phylip'
    else:
        raise ValueError("Invalid input format for {}. Ending determines filetype and should be either .phy or .fa"
                         .format(conf.msa))
    with open(conf.msa) as fh:
        msa = next(Bio.AlignIO.parse(fh, format))
    logger.info("loaded msa with {} sequences of length {}".format(len(msa), msa.get_alignment_length()))

    prepare_r2t(msa, conf.out, conf.min_divergence)

