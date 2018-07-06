import pysam
import numpy as np


class SeqCompleteness(object):

    def __init__(self, mapped_ref=None, tested_ref=None):
        self.seq_completeness = {}

        if mapped_ref:
            self.map_ref_records = self._get_og_dict(mapped_ref)
        else:
            self.map_ref_records = None

        if tested_ref:
            self.ref_records = self._get_og_dict(tested_ref)
        else:
            self.ref_records = None

    def get_seq_completeness(self, records):
        for record in records:
            self.seq_completeness[
                record.id] = self._get_single_seq_completeness(record)

    def _get_single_seq_completeness(self, mapped_record, gene_code='dna'):
        """
        Calculate single sequence completeness using the number of dna or aa
        positions that are not n/X divided by either
        length of sequence or full length or reference
        :param mapped_record: sequence record that was produced by mapping
        :param gene_code: dna or aa
        :return: tuple with partial seq completeness computed using just the
            mapped_record itself and ref_seq_completeness computed
            using also t
        """

        map_ref_record = self.map_ref_records[self._get_og_id(mapped_record.id)]
        map_ref_seq = str(map_ref_record.seq).upper()
        map_seq = str(mapped_record.seq).upper()
        if self.ref_records and self._get_og_id(mapped_record.id) in \
                self.ref_records.keys():
            ref_record = self.ref_records[self._get_og_id(mapped_record.id)]
            ref_seq = str(ref_record.seq).upper()
        else:
            ref_seq = map_ref_seq

        if gene_code is 'dna':
            ref_seq_len = len(ref_seq)
            map_seq_len = len(map_ref_seq)
            non_n_len = len(map_ref_seq) - str(map_seq).count('N')
            map_seq_completeness = non_n_len / map_seq_len
            ref_seq_completeness = non_n_len / ref_seq_len
        elif gene_code is 'aa':
            ref_seq_len = len(ref_seq)
            map_seq_len = len(map_seq)
            non_n_len = len(map_seq) - str(map_seq).count('X')
            map_seq_completeness = non_n_len / map_seq_len
            ref_seq_completeness = non_n_len / ref_seq_len
        return (map_seq_completeness, ref_seq_completeness,
                non_n_len, map_seq_len, ref_seq_len)

    def _get_og_dict(self, ref_og):
        dna_dict = {}
        for record in ref_og:
            if '_' in record.id:
                split_id = record.id.split("_")
                tmp = split_id[0]+"_"+split_id[1]
                record.id = tmp
                og_id = split_id[1]

            dna_dict[og_id] = record
        return dna_dict

    def _get_og_id(self, id):
        split_id = id.split("_")
        # return split_id[0]+"_"+split_id[1]
        return split_id[1]

    def _get_gene_id(self, id):
        split_id = id.split("_")
        return split_id[0]

    def add_seq_completeness(self, ref, seq_completeness):
        self.seq_completeness[ref] = seq_completeness

    def write_seq_completeness(self, file_name):
        out_text = ''
        header = '#species,og,gene_id,map_seq_completeness,' \
                 'ref_seq_completeness,inferred_len,given_len,ref_len\n'
        out_text += header
        for key, value in self.seq_completeness.items():
            species = key[0:5]
            og = key.split("_")[-1]
            gene_id = key.split("_")[0]
            seq_completeness = value
            line = species + "," + og + "," + gene_id + "," + \
                str(seq_completeness[0]) + "," + str(seq_completeness[1]) + \
                "," + str(seq_completeness[2]) + "," + \
                str(seq_completeness[3]) + "," + \
                str(seq_completeness[4]) + "\n"
            out_text += line

        with open(file_name, "w") as myfile:
            myfile.write(out_text)
