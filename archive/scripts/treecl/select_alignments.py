from Bio import AlignIO
import tqdm, os, glob

def run(afolder, ofolder, min_species):
    if not os.path.exists(ofolder):
        os.makedirs(ofolder)
    for f in tqdm.tqdm(glob.glob(os.path.join(afolder, '*.fa'))):
        if os.path.getsize(f) > 0:
            try:
                msa = AlignIO.read(f, "phylip-relaxed")
            except ValueError:
                msa = AlignIO.read(f, "fasta")
        if len(msa) >= min_species:
            align_output = open(os.path.join(ofolder, os.path.basename(f).split(".")[0]+".phy"), "w")
            AlignIO.write(msa, align_output, "phylip-relaxed")
            align_output.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Transform OrthogroupsFasta into marker_genes""")
    parser.add_argument('--afolder', default=None, required=True,
                        help='[Default is none] Folder that contains alignments'
                        'without spaces, e.g. XXX,YYY,AAA.')
    parser.add_argument('--ofolder', default='alignments_selected', required=True,
                        help='[Default is current directory] Path to '
                        'output directory.')
    parser.add_argument('--min_species', type=int, default=0,
                            help='Min number of species in selected '
                                 'alignments. ')

    conf = parser.parse_args()

    run(conf.afolder, conf.ofolder, conf.min_species)
