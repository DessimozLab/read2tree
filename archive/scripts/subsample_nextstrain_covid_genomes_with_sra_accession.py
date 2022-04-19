import itertools
import lzma
import csv
import random


def get_sra_datasets(fn):
    with lzma.open(fn, "rt", newline="") as fh:
        reader = csv.DictReader(fh, dialect="excel-tab")
        for row in reader:
            if row["sra_accession"] not in ('', '?'):
                yield row


def subsample(metafile, nr_per_clade):
    sra = sorted(get_sra_datasets(metafile), key=lambda x: x["Nextstrain_clade"])
    sub = []
    for clade, samples in itertools.groupby(sra, key=lambda x: x["Nextstrain_clade"]):
        if clade == "": 
            continue
        samples = list(samples)
        print(f"{clade}: {len(samples)}")
        sub.extend(random.sample(samples, min(nr_per_clade, len(samples))))
    return sub

def write(outfn, sub):
    with open(outfn,'w') as fout:
        w = csv.DictWriter(fout, fieldnames=sub[0].keys(), dialect="excel-tab")
        w.writeheader()
        w.writerows(sub)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="subsample nextstrain samples from all clades that contain sra accession ids")
    parser.add_argument("--out", required=True, help="path to output file")
    parser.add_argument("--nr-per-clade", default=2, type=int, help="number of samples to use per nextstrain clade. [default: 2]")
    parser.add_argument("metafile", help="metadata.tsv.xz file from nextstrain, e.g. https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.xz")
    conf = parser.parse_args()

    subset = subsample(conf.metafile, conf.nr_per_clade)
    write(conf.out, subset)
