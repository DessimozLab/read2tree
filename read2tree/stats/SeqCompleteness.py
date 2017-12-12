import pysam
import numpy as np


class SeqCompleteness(object):

    def __init__(self, reference=None):
        self.seq_completeness = {}

        if reference:
            self.ref_records = self._get_og_dict(reference)

    def get_seq_completeness(self, records):
        for record in records:
            self.seq_completeness[self._get_clean_id(record.id)] = self._get_single_seq_completeness(record)

    def _get_single_seq_completeness(self, mapped_record, gene_code='dna'):
        """
        Calculate single sequence completeness using the number of dna or aa positions that are not n/X divided by either
        length of sequence or full length or reference
        :param mapped_record: sequence record that was produced by mapping
        :param gene_code: dna or aa
        :return: tuple with partial seq completeness computed using just the mapped_record itself
            and full_seq_completeness computed using also t
        """
        ref_record = self.ref_records[self._get_clean_id(mapped_record.id)]
        if gene_code is 'dna':
            full_seq_len = len(ref_record.seq)
            seq_len = len(mapped_record.seq)
            non_n_len = len(mapped_record.seq) - str(mapped_record.seq).count('n')
            seq_completeness = non_n_len / seq_len
            full_seq_completeness = non_n_len / full_seq_len
        elif gene_code is 'aa':
            full_seq_len = len(ref_record.seq)
            seq_len = len(mapped_record.seq)
            non_n_len = len(mapped_record.seq) - str(mapped_record.seq).count('X')
            seq_completeness = non_n_len / seq_len
            full_seq_completeness = non_n_len / full_seq_len

        return (seq_completeness, full_seq_completeness)

    def _get_og_dict(self, ref_og):
        dna_dict = {}
        for record in ref_og.dna:
            if '_' in record.id:
                split_id = record.id.split("_")
                tmp = split_id[0]+"_"+split_id[1]
                record.id = tmp

            dna_dict[record.id] = record
        return dna_dict

    def _get_clean_id(self, id):
        split_id = id.split("_")
        return split_id[0]+"_"+split_id[1]

    def add_seq_completeness(self, ref, seq_completeness):
        self.seq_completeness[ref] = seq_completeness

    def write_seq_completeness(self, file_name):
        out_text = ''
        header = '#species,og,gene_id,partial_seq_completeness,full_seq_completeness\n'
        out_text += header
        for key, value in self.seq_completeness.items():
            species = key[0:5]
            og = key.split("_")[-1]
            gene_id = key.split("_")[0]
            seq_completeness = value
            line = species + "," + og + "," + gene_id + "," + str(seq_completeness[0]) + "," + str(seq_completeness[1]) + "\n"
            out_text += line

        with open(file_name, "w") as myfile:
            myfile.write(out_text)

    def read_coverage_from_file(self, file_name):
        raise NotImplementedError

    def _get_gene_coverage(self, mybam, ref):
        """

        :param mybam: bam_file object from pysam
        :param ref: the gene_id reference to pileup the the number of reads per column
        :return: average coverage per gene
        """
        column_coverage = []
        for pileupcolumn in mybam.pileup(ref, 0, 100000):
            column_coverage.append(pileupcolumn.n)
        np_column_coverage = np.array(column_coverage)
        return [np.mean(np_column_coverage), np.std(np_column_coverage)]
