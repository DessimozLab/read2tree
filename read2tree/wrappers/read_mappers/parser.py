import logging
import pysam
from pyparsing import Suppress, SkipTo, Word, Regex, Literal, OneOrMore, \
    Group, LineEnd, CharsNotIn, nums, alphanums, ParseException


logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
INT = Word(nums).setParseAction(lambda x: int(x[0]))
WORD = Word(alphanums + '_')
SPACEDWORD = Word(alphanums+' _')


class NGMParser(object):
    """
    Simple prottest result parser.
    """

    def __init__(self):
        self.VAILD_PAIRS = Literal('[MAIN] Valid pairs found:')
        self.INSERT_SIZE = Literal('[MAIN] Estimated insert size:')
        self.COMPUTED_ALIGNMENTS = Literal('[MAIN] Alignments computed:')
        # These are all the models that are possible to be tested using phyml
        self.vp = Suppress(SkipTo(self.VAILD_PAIRS)) + Suppress(self.VAILD_PAIRS) + WORD
        self.inssize = Suppress(SkipTo(self.INSERT_SIZE)) + Suppress(self.INSERT_SIZE) + FLOAT
        self.compalign = Suppress(SkipTo(self.COMPUTED_ALIGNMENTS)) + Suppress(self.COMPUTED_ALIGNMENTS) + FLOAT

    def parse(self, stdout):
        try:
            valid_pairs = self.vp.parseString(stdout).asList()[0]
            insert_size = self.inssize.parseString(stdout).asList()[0]
            computed_alignments = self.compalign.parseString(stdout).asList()[0]

        except ParseException as err:
            logger.error(err)

        return valid_pairs, insert_size, computed_alignments

    def to_dict(self, file, stdout):
        valid_pairs, insert_size, computed_alignments = self.parse(stdout)
        samfile = pysam.AlignmentFile(file, "r")
        result = {'file': file,
                  'valid_pairs': valid_pairs,
                  'insert_size': insert_size,
                  'computed_alignments': computed_alignments,
                  'sam': samfile}

        return result

class NGMLRParser(object):
    """
    Simple prottest result parser.
    for the following example output line:
        Processed: 75400 (0.00), R/S: 60.15, RL: 7675, Time: 3.00 11.00 10.07, Align: 1.00, 310, 3.04
        Done (77 reads mapped (0.10%), 75323 reads not mapped, 75402 lines written)(elapsed: 20m, 0 r/s)
    """

    def __init__(self):
        self.TOTAL_MAPPED_READS = Literal('Done (')
        self.TOTAL_READS = Literal('Processed: ')
        # These are all the models that are possible to be tested using phyml
        self.tmr = Suppress(SkipTo(self.TOTAL_MAPPED_READS)) + Suppress(self.TOTAL_MAPPED_READS) + FLOAT
        self.tr = Suppress(SkipTo(self.TOTAL_READS)) + Suppress(self.TOTAL_READS) + FLOAT

    def parse(self, stdout):
        try:
            total_mapped_reads = self.tmr.parseString(stdout).asList()[0]
            total_reads = self.tr.parseString(stdout).asList()[0]

        except ParseException as err:
            logger.error(err)

        return total_mapped_reads, total_reads

    def to_dict(self, file, stdout):
        total_mapped_reads, total_reads = self.parse(stdout)
        samfile = pysam.AlignmentFile(file, "r")
        result = {'file': file,
                  'total_mapped_reads': total_mapped_reads,
                  'total_reads': total_reads,
                  'sam': samfile}

        return result






