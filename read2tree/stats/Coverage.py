import pysam

class Coverage(object):

    def __init__(self):
        self.coverage = {}

    def get_coverage_bam(self, file_name):
        mybam = pysam.AlignmentFile(file_name, 'rb')
        for ref in mybam.references:
            self.coverage[ref] = self._get_gene_coverage(mybam, ref)
        #return self.coverages

    def write_coverage_bam(self, file_name):
        out_text = ''
        header = '#species,og,gene_id,coverage\n'
        out_text += header
        for key, value in self.coverage.items():
            species = key[0:5]
            og = key.split("_")[-1]
            gene_id = key.split("_")[0]
            coverage = value
            line = species + "," + og + "," + gene_id + "," + str(coverage) + "\n"
            out_text += line


        with open(file_name, "w") as myfile:
            myfile.write(out_text)

    def read_coverage_from_file(self, file_name):
        raise NotImplementedError


    def _get_gene_coverage(self, mybam, ref):
        absolute_coverage = 0
        num_positions = 0
        relative_coverage = 0
        for pileupcolumn in mybam.pileup(ref, 0, 100000):
            absolute_coverage += pileupcolumn.n
            num_positions += 1
        if num_positions != 0:
            relative_coverage = absolute_coverage/num_positions

        return relative_coverage