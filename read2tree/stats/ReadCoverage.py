

class Coverage(object):

    def __init__(self, file_name):
        self.coverage = self._get_coverage_bam(file_name)

    def _get_coverage_bam(self, file_name):
        mybam = pysam.AlignmentFile(file_name, 'rb')
        for ref in mybam.references:
            self.coverage[ref] = self._get_gene_coverage(mybam, ref)


    def _get_gene_coverage(self, mybam, ref):
        coverage = 0
        num_positions = 0
        relative_coverage = 0
        for pileupcolumn in mybam.pileup(ref, 0, 100000):
            absolute_coverage += pileupcolumn.n
            num_positions += 1
        if num_positions != 0:
            relative_coverage = absolute_coverage/num_positions

        return relative_coverage

    def _write_coverage_bam(self):
        raise NotImplementedError


class GeneCoverage(object):

    def __init__(self):
