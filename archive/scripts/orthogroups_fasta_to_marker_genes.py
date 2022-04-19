from Bio.SeqIO.FastaIO import FastaWriter
from Bio import SeqIO
import tqdm, os, glob, re
from xml.dom import minidom



def _find_index_substring(ids, search_string, tmp_list):
    best_index = None
    max_occurence = 0
    tmp_ids = [re.sub(r'\..*', '', tmp) for tmp in tmp_list]
    use_ids = [re.sub(r'\W+', '', tmp_id) for tmp_id in tmp_ids]
    index = [i for i, s in enumerate(ids) if search_string in s]
    for i in index:
        string_occurence = len([k for k in use_ids if k in ids[i]])
        if string_occurence > max_occurence:
            best_index = i
            max_occurence = string_occurence
    if best_index:
        return best_index
    else:
        return None


def _get_all_ids(f_orthoxml):
    all_prot_ids = []
    xmldoc = minidom.parse(f_orthoxml)
    itemlist = xmldoc.getElementsByTagName('gene')
    print(" --- loading all protids ---")
    for s in tqdm.tqdm(itemlist):
        tmp = s.attributes['protId'].value
        all_prot_ids.append(tmp)
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


def _get_species_id(record):
    if '[' in record.description and ']' in record.description:
        return record.description[record.description.find(
            "[")+1:record.description.find("]")]
    else:
        return record.id[0:5]

def run(orthogroups_fasta_folder, orthogroups_xml, output_path, min_species):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    all_prot_ids = _get_all_ids(orthogroups_xml)
    for f in tqdm.tqdm(glob.glob(os.path.join(orthogroups_fasta_folder, '*.fa'))):
        records = list(SeqIO.parse(f, 'fasta'))
        if len(records) >= min_species:
            for rec in records:
                sp_id = _get_species_id(rec)
                tmp_lst = rec.description.split()
                if sp_id not in tmp_lst[0]:
                    tmp = tmp_lst[-2]
                    tmp_id = re.sub(r'\..*', '', tmp)
                    use_id = re.sub(r'\W+', '', tmp_id)
                    new_id = _find_index_substring(all_prot_ids, use_id, tmp_lst)
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
    parser.add_argument('--min_species', type=int, default=None,
                            help='Min number of species in selected '
                                 'orthologous groups. If not selected it will be '
                                 'estimated such that around 1000 OGs '
                                 'are available.')

    conf = parser.parse_args()

    run(conf.ofasta, conf.oxml, conf.ofolder, conf.min_species)
