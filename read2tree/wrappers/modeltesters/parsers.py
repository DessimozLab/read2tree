import logging
import dendropy as dpy
from pyparsing import Suppress, SkipTo, Word, Regex, Literal, OneOrMore, Group, LineEnd, CharsNotIn, nums, alphanums, \
    ParseException

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
INT = Word(nums).setParseAction(lambda x: int(x[0]))
WORD = Word(alphanums + '_')
SPACEDWORD = Word(alphanums+' _')


class ProtTestParser(object):
    """
    Simple prottest result parser.
    """

    def __init__(self):
        self.MODEL = Regex(r'Best model according to\s+')
        # These are all the models that are possible to be tested using phyml
        self.model = OneOrMore(Group(Suppress(SkipTo(self.MODEL)) + Suppress(self.MODEL) + WORD + Suppress(":") + WORD))

    def parse(self, s):
        model = None
        try:
            model = self.model.parseString(s).asList()
        except ParseException as err:
            logger.error(err)

        return model

    def to_dict(self, stats_filename):
        result = {}
        model = self.parse(stats_filename)
        try:
            for mg in model:
                result[mg[0]] = mg[1]
        except IOError as err:
            logger.error(err)
            return

        return result




