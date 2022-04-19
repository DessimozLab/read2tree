from Bio.SeqIO.FastaIO import FastaWriter
from Bio import SeqIO
import tqdm, os, glob

def _oma_replace(row):
    if 'OMA0000' in row:
        return 'OMA0000'
    elif 'OMA000' in row:
        return 'OMA000'
    elif 'OMA00' in row:
        return 'OMA00'
    elif 'OMA0' in row:
        return 'OMA0'
    elif 'OMA' in row:
        return 'OMA'


def _get_all_ids(orthogroups_txt):
    with open(orthogroups_txt) as f:
        lines = f.readlines()
    x = []
    for l in lines:
        if '#' not in l:
            x.append(l.rstrip("\n").split("\t"))
    og_dic = {}
    for r in x:
        tmp = r[0].replace(_oma_replace(r[0]), 'OG')
        r[0] = tmp
        og_dic[tmp] = {i[0:5]: i[6:] for i in r[1:]}
    return og_dic


def _write(file, value):
    """
    Write output to fasta file
            :param file: file and location of outputfile
            :param value:
            :return:
    """
    handle = open(file, "w")
    writer = FastaWriter(handle, wrap=None)
    writer.write_file(value)
    handle.close()


def _get_species_id(record):
    if '[' in record.description and ']' in record.description:
        return record.description[record.description.find(
            "[")+1:record.description.find("]")]
    else:
        return record.id[0:5]

def run(orthogroups_fasta_folder, og_dic, output_path, min_species):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for f in tqdm.tqdm(glob.glob(os.path.join(orthogroups_fasta_folder, '*.fa'))):
        new_name_dic = og_dic[os.path.basename(f).split(".")[0]]
        records = list(SeqIO.parse(f, 'fasta'))
        if len(records) >= min_species:
            for rec in records:
                sp_id = _get_species_id(rec)
                new_id = new_name_dic[sp_id].split()[0]
                rec.id = new_id
                rec.description = new_name_dic[sp_id].replace(new_id, "") + " [" + sp_id + "]"
            output_file = os.path.join(output_path,
                                           os.path.basename(f))
            _write(output_file, records)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Transform OrthogroupsFasta into marker_genes""")
    parser.add_argument('--ogroups', default=None,
                        help='[Default is none] Remove species present '
                        'in data set after mapping step completed to '
                        'build OGs. Input is comma separated list '
                        'without spaces, e.g. XXX,YYY,AAA.')
    parser.add_argument('--ofolder', default='marker_genes', required=True,
                        help='[Default is current directory] Path to '
                        'output directory.')
    parser.add_argument('--ofasta', default='.', required=True,
                        help='[Default is current directory] Path to '
                        'output directory.')
    parser.add_argument('--min_species', type=int, default=None,
                            help='Min number of species in selected '
                                 'orthologous groups. If not selected it will be '
                                 'estimated such that around 1000 OGs '
                                 'are available.')

    conf = parser.parse_args()
    og_dic = _get_all_ids(conf.ogroups)

    run(conf.ofasta, og_dic, conf.ofolder, conf.min_species)
