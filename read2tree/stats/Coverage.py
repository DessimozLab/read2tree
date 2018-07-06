import pysam
import numpy as np


class Coverage(object):

    def __init__(self):
        self.coverage = {}

    def get_coverage_bam(self, file_name):
        mybam = pysam.AlignmentFile(file_name, 'rb')
        for ref in mybam.references:
            self.coverage[self._get_clean_id(ref)] \
                = self._get_gene_coverage(mybam, ref)

    def _get_clean_id(self, id):
        id = id.split(" ")[0]
        id = id.split("_")
        return id[0]+"_"+id[1]

    def add_coverage(self, ref, coverage):
        self.coverage[ref] = coverage

    def write_coverage_bam(self, file_name):
        out_text = ''
        header = '#species,og,gene_id,coverage,std\n'
        out_text += header
        for key, value in self.coverage.items():
            species = key[0:5]
            og = key.split("_")[-1]
            gene_id = key.split("_")[0]
            coverage = value
            line = species + "," + og + "," + gene_id + "," + \
                str(coverage[0]) + "," + str(coverage[1]) + "\n"
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
