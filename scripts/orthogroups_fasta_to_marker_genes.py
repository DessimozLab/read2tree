from Bio.SeqIO.FastaIO import FastaWriter
from Bio import SeqIO
import tqdm, os, glob
from xml.dom import minidom



def _find_index_substring(ids, search_string):
    index = [i for i, s in enumerate(ids) if search_string in s]
    if index:
        return index[0]
    else:
        return None


def _get_all_ids(f_orthoxml):
    all_prot_ids = []
    xmldoc = minidom.parse(f_orthoxml)
    itemlist = xmldoc.getElementsByTagName('gene')
    print(" --- loading all protids ---")
    for s in tqdm.tqdm(itemlist):
        all_prot_ids.append(s.attributes['protId'].value)
    return all_prot_ids


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


def run(orthogroups_fasta_folder, orthogroups_xml, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    all_prot_ids = _get_all_ids(orthogroups_xml)
    for f in tqdm.tqdm(glob.glob(os.path.join(orthogroups_fasta_folder, '*.fa'))):
        records = list(SeqIO.parse(f, 'fasta'))
        for rec in records:
            new_id = _find_index_substring(all_prot_ids, rec.id)
            if new_id:
                rec.id = all_prot_ids[new_id]
                new_description = rec.description.split()[-1]
                rec.description = new_description
                rec.name = ''
        output_file = os.path.join(output_path,
                                   os.path.basename(f))
        _write(output_file, records)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Transform OrthogroupsFasta into marker_genes""")
    parser.add_argument('--oxml', default=None,
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

    conf = parser.parse_args()

    run(conf.ofasta, conf.oxml, conf.ofolder)
