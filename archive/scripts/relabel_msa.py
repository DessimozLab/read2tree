import Bio.AlignIO
import csv


def load_oma_species(fn):
    with open(fn, 'rt') as fh:
        reader = csv.reader((l for l in fh if not l.startswith('#')), dialect="excel-tab")
        mapping = {row[0]: row[2].replace(' ','_') + "__" + row[1] for row in reader}
    return mapping


def load_nextstrain_metadata(fn):
    with open(fn, 'rt') as fh:
        reader = csv.DictReader(fh, dialect="excel-tab")
        mapping = {row['sra_accession']: row['sra_accession'] + "__" + row['strain'].replace(' ','_') + "__" + row['Nextstrain_clade'].replace(' ','_').replace('(','[').replace(')',']') + row['date']
                   for row in reader}
    return mapping


def update_msa_ids(msa_path, new_path, mapping, format="phylip-relaxed"):
    msa = Bio.AlignIO.read(msa_path, format=format)
    for rec in msa:
        rec.id = mapping.get(rec.id, rec.id)
    Bio.AlignIO.write(msa, new_path, format=format)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="update labels of sequence ids")
    parser.add_argument('--oma-map', help="path to the oma-species.txt file to change 5letter codes with scientific names")
    parser.add_argument('--nextstrain', help="path to the nextstrain metadata file with the sra accessions")
    parser.add_argument('--msa-format', help="format of the msa. if not set, it will be guessed based on file extension")
    parser.add_argument('--out', required=True, help="Path to the output filename")
    parser.add_argument('msa', help="Path to the input msa filename")

    conf = parser.parse_args()
    mapping = {}
    if conf.oma_map:
        mapping.update(load_oma_species(conf.oma_map))
    if conf.nextstrain:
        mapping.update(load_nextstrain_metadata(conf.nextstrain))

    if conf.msa_format is None:
        conf.msa_format = "phylip-relaxed" if conf.msa.endswith('.phy') else "fasta"
    update_msa_ids(conf.msa, conf.out, mapping, format=conf.msa_format)

