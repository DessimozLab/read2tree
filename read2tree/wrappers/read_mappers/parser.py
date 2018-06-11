import logging
import pysam
from pyparsing import Suppress, SkipTo, Word, Regex, Literal, OneOrMore, \
    Group, LineEnd, CharsNotIn, nums, alphanums, ParseException


logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
INT = Word(nums).setParseAction(lambda x: int(x[0]))
WORD = Word(alphanums + '_-%. ')
SPACEDWORD = Word(alphanums+' _')


class NGMParser(object):
    """
    Simple prottest result parser.
    [MAIN] Done (15778 reads mapped (4.14%), 365184 reads not mapped, 15778 lines written)(elapsed: 73.131973s)
    """

    def __init__(self):
        self.READS_MAPPED = Literal('[MAIN] Done (')
        self.TOTAL_READS = Regex(r'\[MAIN\] Done \(\d+ reads mapped \(\d+\.\d+\%\), ')
        self.MAPPING_TIME = Literal('elapsed: ')
        self.rm = Suppress(SkipTo(self.READS_MAPPED)) + Suppress(self.READS_MAPPED) + INT
        self.tr = Suppress(SkipTo(self.TOTAL_READS)) + Suppress(self.TOTAL_READS) + INT
        self.mt = Suppress(SkipTo(self.MAPPING_TIME)) + Suppress(self.MAPPING_TIME) + FLOAT

    def parse(self, stdout):
        try:
            reads_mapped = self.rm.parseString(stdout).asList()[0]
            total_reads = self.tr.parseString(stdout).asList()[0]
            mapping_time = self.mt.parseString(stdout).asList()[0]
        except ParseException as err:
            print(stdout)
            logger.error(err)
        else:
            return reads_mapped, total_reads, mapping_time

    def to_dict(self, file, stdout):
        try:
            reads_mapped, total_reads, mapping_time = self.parse(stdout)
        except UnboundLocalError:
            reads_mapped = None
            total_reads = None
            mapping_time = None
            pass
        samfile = pysam.AlignmentFile(file, "r")
        result = {'file': file,
                  'reads_mapped': reads_mapped,
                  'total_reads': total_reads,
                  'mapping_time': mapping_time,
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
        self.tmr = Suppress(SkipTo(self.TOTAL_MAPPED_READS)) + \
            Suppress(self.TOTAL_MAPPED_READS) + FLOAT
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
                  'reads_mapped': total_mapped_reads,
                  'total_reads': total_reads,
                  'sam': samfile}

        return result
