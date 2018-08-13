from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped


def get_alignment(file, species_list):
    keep_species = []
    alignment = AlignIO.read(file, 'phylip-relaxed')
    for i, record in enumerate(alignment):
        if record.id not in species_list:
            keep_species.append(record)
    return MultipleSeqAlignment(keep_species, Gapped(IUPAC.protein, "-"))


def write_alignment(output, alignment):
    AlignIO.write(alignment, output, 'phylip-relaxed')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Remove species from given alignment.""")
    parser.add_argument('-s', '--remove_species', default=None,
                        help='[Default is none] Remove species present '
                        'in data set after mapping step completed to '
                        'build OGs. Input is comma separated list '
                        'without spaces, e.g. XXX,YYY,AAA.')
    parser.add_argument('-o', '--output', default='.', required=True,
                        help='[Default is current directory] Path to '
                        'output directory.')
    parser.add_argument('-i', '--input', default='.', required=True,
                        help='[Default is current directory] Path to '
                        'output directory.')

    conf = parser.parse_args()

    new_alignment = get_alignment(conf.input, conf.remove_species)
    write_alignment(conf.output, new_alignment)
