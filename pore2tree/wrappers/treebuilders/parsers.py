import logging
import dendropy as dpy
from pyparsing import Suppress, SkipTo, Word, Regex, Literal, OneOrMore, \
    Group, LineEnd, CharsNotIn, nums, alphanums, ParseException

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
INT = Word(nums).setParseAction(lambda x: int(x[0]))
WORD = Word(alphanums + '_')
SPACEDWORD = Word(alphanums + ' _')


class PhymlParser(object):
    """
    Simple phyml result parser. Assumes one of the standard models
    for nucleotide analyses.
    """

    def __init__(self):
        self.MODEL_LABEL = Regex(r'Model of.*substitution:\s+')
        self.ALPHA_LABEL = Regex(r'Gamma shape parameter:\s+')
        self.LNL_LABEL = Regex(r'Log-likelihood:\s+')
        self.F_LABEL = Regex(r'f\(([ACGT])\)=\s+')
        self.R_LABEL = Regex(r'[ACGT]\s+<->\s+[ACGT]\s+')
        self.TSTV_LABEL = Regex(r'Transition/transversion ratio.*:\s+')
        self.model = Suppress(SkipTo(self.MODEL_LABEL)) + \
                     Suppress(self.MODEL_LABEL) + WORD
        self.lnl = Suppress(SkipTo(self.LNL_LABEL)) + \
                   Suppress(self.LNL_LABEL) + FLOAT
        self.alpha = Suppress(SkipTo(self.ALPHA_LABEL)) + \
                     Suppress(self.ALPHA_LABEL) + FLOAT
        self.common = self.model + self.lnl + self.alpha
        self.tstv = OneOrMore(Suppress(SkipTo(self.TSTV_LABEL)) +
                              Suppress(self.TSTV_LABEL) + FLOAT)
        self.freq = OneOrMore(
            Suppress(SkipTo(self.F_LABEL)) + Suppress(self.F_LABEL) + FLOAT)
        self.rates = OneOrMore(
            Suppress(SkipTo(self.R_LABEL)) + Suppress(self.R_LABEL) + FLOAT)
        self.gtr_specific = Group(self.freq) + Group(self.rates)
        self.hky_specific = Group(self.tstv) + Group(self.freq)

    def parse(self, filename):
        model = None
        alpha = None
        lnl = None
        freq = None
        rates = None

        with open(filename) as fl:
            s = fl.read()

        try:
            model, lnl, alpha = self.common.parseString(s).asList()

        except ParseException as err:
            logger.error(err)

        if model == 'JC69':
            freq = [0.25, 0.25, 0.25, 0.25]
            rates = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        elif model == 'K80':
            freq = [0.25, 0.25, 0.25, 0.25]
            try:
                tstv = self.tstv.parseString(s).asList()
            except ParseException as err:
                logger.error(err)

            rates = [1.0, tstv[0], 1.0, 1.0, tstv[0], 1.0]

        elif model == 'F81':
            try:
                freq = self.freq.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
            rates = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        elif model == 'F84' or model == 'HKY85' or model == 'TN93':
            parser = Group(self.tstv) + Group(self.freq)
            try:
                tstv, freq = parser.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
            if model == 'TN93':
                rates = [1.0, tstv[0], 1.0, 1.0, tstv[1], 1.0]
            else:
                rates = [1.0, tstv[0], 1.0, 1.0, tstv[0], 1.0]

        elif model == 'GTR':
            parser = Group(self.freq) + Group(self.rates)
            try:
                freq, rates = parser.parseString(s).asList()
            except ParseException as err:
                logger.error(err)

        return model, alpha, lnl, freq, rates

    def to_dict(self, stats_filename, tree_filename):
        model, alpha, lnl, freq, rates = self.parse(stats_filename)
        try:
            with open(tree_filename) as treefl:
                nwk_tree = treefl.read().rstrip()
                tree = dpy.Tree.get(data=nwk_tree, schema='newick')
        except IOError as err:
            logger.error(err)
            return

        result = {'likelihood': lnl,
                  'alpha': alpha,
                  'frequencies': freq,
                  'rates': rates,
                  'model': model,
                  'tree': tree}
        return result


class RaxmlParser(object):
    def __init__(self):
        self.ALPHA_LABEL = Regex(r'alpha\[\d+\]:')
        self.LNL_LABEL = Literal('Final GAMMA-based Score of best tree')
        self.FRQ_LABEL = Regex(r'Base frequencies: (?=\d+)') ^ \
                         Regex(r'ML estimate base freqs\[\d+\]:')
        self.NAMES_LABEL = Regex(r'Partition: \d+ with name:\s+')
        self.RATES_LABEL = Regex(r'rates\[\d+\].+?:')
        self.MODEL_LABEL = Literal('Substitution Matrix:')
        self.alpha = OneOrMore(Suppress(SkipTo(self.ALPHA_LABEL)) +
                               Suppress(self.ALPHA_LABEL) + FLOAT)
        self.lnl = Suppress(SkipTo(self.LNL_LABEL)) + \
                   Suppress(self.LNL_LABEL) + FLOAT
        self.frq = OneOrMore(Group(Suppress(SkipTo(self.FRQ_LABEL)) +
                                   Suppress(self.FRQ_LABEL) + OneOrMore(FLOAT)))
        self.names = OneOrMore(
            Suppress(SkipTo(self.NAMES_LABEL)) + Suppress(self.NAMES_LABEL) +
            CharsNotIn('\n') + Suppress(LineEnd()))
        self.rates = OneOrMore(
            Group(Suppress(SkipTo(self.RATES_LABEL)) +
                  Suppress(self.RATES_LABEL) + OneOrMore(FLOAT)))
        self.model = Suppress(SkipTo(self.MODEL_LABEL)) + \
                     Suppress(self.MODEL_LABEL) + WORD

        MODEL_LABEL = Literal('Substitution Matrix:')
        SCORE_LABEL = Literal('Final GAMMA  likelihood:')
        BOOT_SCORE_LABEL = Literal('Final ML Optimization Likelihood:')
        DESC_LABEL = Literal('Model Parameters of Partition')
        NAME_LEADIN = Literal(', Name:')
        DATATYPE_LEADIN = Literal(', Type of Data:')
        ALPHA_LEADIN = Literal('alpha:')
        TREELENGTH_LEADIN = Literal('Tree-Length:')
        RATES_LABEL = Regex(r'rate \w <-> \w:')
        FREQS_LABEL = Regex(r'freq pi\(\w\):')

        likelihood = Suppress(SkipTo(SCORE_LABEL)) + Suppress(
            SCORE_LABEL) + FLOAT
        boot_likelihood = Suppress(SkipTo(BOOT_SCORE_LABEL)) + Suppress(
            BOOT_SCORE_LABEL) + FLOAT
        description = Suppress(SkipTo(DESC_LABEL)) + Suppress(
            DESC_LABEL) + INT + Suppress(
            NAME_LEADIN) + SPACEDWORD + Suppress(DATATYPE_LEADIN) + WORD
        treelen = Suppress(SkipTo(TREELENGTH_LEADIN)) + Suppress(
            TREELENGTH_LEADIN) + FLOAT

        alpha = Suppress(SkipTo(ALPHA_LEADIN)) + Suppress(ALPHA_LEADIN) + FLOAT

        rates = OneOrMore(
            Group(Suppress(SkipTo(RATES_LABEL)) + Suppress(
                RATES_LABEL) + OneOrMore(FLOAT)))
        freqs = OneOrMore(
            Group(Suppress(SkipTo(FREQS_LABEL)) + Suppress(
                FREQS_LABEL) + OneOrMore(FLOAT)))

        # output of running different set of raxml analysis
        self.TC_STOCHBI_LABEL = Literal(
            'Tree certainty under stochastic bipartition '
            'adjustment for this tree:')
        self.RTC_STOCHBI_LABEL = Literal(
            'Relative tree certainty under stochastic bipartition adjustment for this tree:')
        self.TCA_STOCHBI_LABEL = Literal(
            'Tree certainty including all conflicting bipartitions (TCA) under '
            'stochastic bipartition adjustment for this tree:')
        self.RTCA_STOCHBI_LABEL = Literal(
            'Relative tree certainty including all conflicting bipartitions (TCA) '
            'under stochastic bipartition adjustment for this tree:')
        self.TC_UNIBI_LABEL = Literal(
            'Tree certainty under uniform bipartition '
            'adjustment for this tree:')
        self.RTC_UNIBI_LABEL = Literal('Relative tree certainty under uniform '
                                       'bipartition adjustment for this tree:')
        self.TCA_UNIBI_LABEL = Literal(
            'Tree certainty including all conflicting bipartitions (TCA) under '
            'uniform bipartition adjustment for this tree:')
        self.RTCA_UNIBI_LABEL = Literal(
            'Relative tree certainty including all conflicting bipartitions (TCA) '
            'under uniform bipartition adjustment for this tree:')

        self.tc_stochbi = Suppress(SkipTo(self.TC_STOCHBI_LABEL)) + Suppress(
            self.TC_STOCHBI_LABEL) + FLOAT
        self.rtc_stochbi = Suppress(SkipTo(self.RTC_STOCHBI_LABEL)) + Suppress(
            self.RTC_STOCHBI_LABEL) + FLOAT
        self.tca_stochbi = Suppress(SkipTo(self.TCA_STOCHBI_LABEL)) + Suppress(
            self.TCA_STOCHBI_LABEL) + FLOAT
        self.rtca_stochbi = Suppress(
            SkipTo(self.RTCA_STOCHBI_LABEL)) + Suppress(
            self.RTCA_STOCHBI_LABEL) + FLOAT
        self.tc_unibi = Suppress(SkipTo(self.TC_UNIBI_LABEL)) + Suppress(
            self.TC_UNIBI_LABEL) + FLOAT
        self.rtc_unibi = Suppress(SkipTo(self.RTC_UNIBI_LABEL)) + Suppress(
            self.RTC_UNIBI_LABEL) + FLOAT
        self.tca_unibi = Suppress(SkipTo(self.TCA_UNIBI_LABEL)) + Suppress(
            self.TCA_UNIBI_LABEL) + FLOAT
        self.rtca_unibi = Suppress(SkipTo(self.RTCA_UNIBI_LABEL)) + Suppress(
            self.RTCA_UNIBI_LABEL) + FLOAT

        # Use these for flag 'a' option
        self.boot_likelihood = boot_likelihood
        self.freqs = freqs
        self.rates = rates
        self.alpha = alpha
        self.name = description
        self.treelen = treelen

        self._dash_f_e_parser = (Group(OneOrMore(self.model)) +
                                 likelihood +
                                 Group(OneOrMore(Group(
                                     description + alpha + Suppress(
                                         TREELENGTH_LEADIN) +
                                     Suppress(FLOAT) + Group(
                                         OneOrMore(rates)) + Group(
                                         OneOrMore(freqs))))))

    def parse(self, filename, f_flag=None):
        with open(filename) as fl:
            s = fl.read()
        if f_flag is None:
            try:
                alphas = self.alpha.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
                alphas = [None]
            try:
                freqs = self.frq.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
                freqs = [None]
            try:
                names = self.names.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
                names = [None]
            try:
                rates = self.rates.parseString(s).asList()
            except ParseException:
                rates = None
            try:
                lnl = self.lnl.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
                lnl = [0]

            return alphas, freqs, names, rates, lnl
        elif f_flag is 'i':
            try:
                tc_stochbi = self.tc_stochbi.parseString(s).asList()[0]
            except ParseException as err:
                logger.error(err)
                tc_stochbi = [None]
            try:
                rtc_stochbi = self.rtc_stochbi.parseString(s).asList()[0]
            except ParseException as err:
                logger.error(err)
                rtc_stochbi = [None]
            try:
                tca_stochbi = self.tca_stochbi.parseString(s).asList()[0]
            except ParseException as err:
                logger.error(err)
                tca_stochbi = [None]
            try:
                rtca_stochbi = self.rtca_stochbi.parseString(s).asList()[0]
            except ParseException as err:
                logger.error(err)
                rtca_stochbi = [None]
            try:
                tc_unibi = self.tc_unibi.parseString(s).asList()[0]
            except ParseException as err:
                logger.error(err)
                tc_unibi = [None]
            try:
                rtc_unibi = self.rtc_unibi.parseString(s).asList()[0]
            except ParseException as err:
                logger.error(err)
                rtc_unibi = [None]
            try:
                tca_unibi = self.tca_unibi.parseString(s).asList()[0]
            except ParseException as err:
                logger.error(err)
                tca_unibi = [None]
            try:
                rtca_unibi = self.rtca_unibi.parseString(s).asList()[0]
            except ParseException as err:
                logger.error(err)
                rtca_unibi = [None]

            return tc_unibi, rtc_unibi, tca_unibi, rtca_unibi, tc_stochbi, \
                   rtc_stochbi, tca_stochbi, rtca_stochbi

        elif f_flag is 'a':
            try:
                model = self.model.parseString(s).asList()
            except ParseException:
                model = None
            try:
                alphas = self.alpha.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
                alphas = [None]
            try:
                freqs = self.freqs.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
                freqs = [None]
            try:
                name = self.name.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
                name = [None]
            try:
                rates = self.rates.parseString(s).asList()
            except ParseException:
                rates = None
            try:
                likelihood = self.boot_likelihood.parseString(s).asList()
            except ParseException:
                likelihood = None
            try:
                treelength = self.treelen.parseString(s).asList()
            except ParseException:
                treelength = None
            return model, treelength, alphas, freqs, name, rates, likelihood

    def _dash_f_e_to_dict(self, info_filename, tree_filename):
        """
        Raxml provides an option to fit model params to a tree,
        selected with -f e.
        The output is different and needs a different parser.
        """
        with open(info_filename) as fl:
            models, likelihood, partition_params = self._dash_f_e_parser.parseFile(
                fl).asList()

        with open(tree_filename) as fl:
            nwk_tree = fl.read().rstrip()
            tree = dpy.Tree.get(data=nwk_tree, schema='newick')

        d = {'likelihood': likelihood, 'tree': tree, 'partitions': {}}

        for model, params in zip(models, partition_params):
            subdict = {}
            index, name, _, alpha, rates, freqs = params
            subdict['alpha'] = alpha
            subdict['name'] = name
            subdict['rates'] = rates
            subdict['frequencies'] = freqs
            subdict['model'] = model
            d['partitions'][index] = subdict

        return d

    def _dash_f_i_to_dict(self, info_filename, tree_filename):
        """
        Raxml provides an option to fit model params to a tree,
        selected with -f i.
        The output is different and needs a different parser.
        """
        tc_unibi, rtc_unibi, tca_unibi, rtca_unibi, tc_stochbi, rtc_stochbi, tca_stochbi, rtca_stochbi = self.parse(
            info_filename, 'i')

        with open(tree_filename) as fl:
            nwk_tree = fl.read().rstrip()
            tree = dpy.Tree.get(data=nwk_tree, schema='newick')

        d = {'tree': tree, 'tc_unibi': tc_unibi, 'rtc_unibi': rtc_unibi,
             'tca_unibi': tca_unibi,
             'rtca_unibi': rtca_unibi, 'tc_stochbi': tc_stochbi,
             'rtc_stochbi': rtc_stochbi, 'tca_stochbi': tca_stochbi,
             'rtca_stochbi': rtca_stochbi}

        return d

    def _dash_f_a_to_dict(self, info_filename, tree_filename):
        """
        Raxml provides an option to fit model params to a tree,
        selected with -f a.
        The output is different and needs a different parser.
        """
        model, treelength, alphas, freqs, name, rates, likelihood = self.parse(
            info_filename, f_flag='a')

        bs_tree_filename = tree_filename.replace('bestTree',
                                                 'bipartitionsBranchLabels')
        try:
            with open(bs_tree_filename) as fl:
                nwk_tree = fl.read().rstrip()
                tree = dpy.Tree.get(data=nwk_tree, schema='newick')
        except IOError as err:
            logger.error('No tree file - raxml analysis failed')
            return

        result = {'likelihood': likelihood,
                  'model': model,
                  'treelength': treelength,
                  'part_name': name[1],
                  'alpha': alphas,
                  'rates': rates,
                  'frequences': freqs,
                  'bs_tree': nwk_tree,
                  'tree': tree}

        return result

    def to_dict(self, info_filename, tree_filename, dash_f=None):
        """
        Parse raxml output and return a dict
        Option dash_f_e=True will parse the output of a raxml -f e run,
        which has different output
        """
        if dash_f is 'e':
            return self._dash_f_e_to_dict(info_filename, tree_filename)
        elif dash_f is 'i':
            return self._dash_f_i_to_dict(info_filename, tree_filename)
        elif dash_f is 'a':
            return self._dash_f_a_to_dict(info_filename, tree_filename)
        else:
            return self._to_dict(info_filename, tree_filename)

    def _to_dict(self, info_filename, tree_filename):
        alpha, freqs, names, rates, lnl = self.parse(info_filename)
        try:
            with open(tree_filename) as fl:
                nwk_tree = fl.read().rstrip()
                tree = dpy.Tree.get(data=nwk_tree, schema='newick')
        except IOError as err:
            logger.error('No tree file - raxml analysis failed')
            return
        n_parts = len(alpha)
        assert len(freqs) == n_parts
        assert len(names) == n_parts
        if rates is not None:
            assert len(rates) == n_parts

        result = {'likelihood': lnl[0],
                  'partitions': {},
                  'tree': tree}

        for i in range(n_parts):
            subdict = {'alpha': alpha[i], 'frequencies': freqs[i],
                       'name': names[i]}
            if rates is not None:
                subdict['rates'] = rates[i]
            result['partitions'][i] = subdict

        return result


class IqtreeParser(object):
    def __init__(self):
        """
        Put stuff here (regex etc.)
        """

    def parse(self, filename):
        """
        Extract based on strings defined in __init__
        """

    def to_dict(self, tree_filename):

        """
        Add stuff for results dictionary here
        """

        try:
            with open(tree_filename) as fl:
                tree = fl.read()
        except IOError as err:
            logger.error(err)
            return

        tree = tree.rstrip()
        result = {'tree': tree}

        return result


class FasttreeParser(object):
    def __init__(self):
        """
        Put stuff here (regex etc.)
        """
        self.tree = None
        self.likelihood = None

    def parse(self, tree, other):
        """
        Extract based on strings defined in __init__
        """
        self.tree = dpy.Tree.get(data=tree.rstrip(), schema='newick')
        for line in other.split('\n'):
            if 'Optimize all lengths:' in line:
                self.likelihood = float(line.split()[-3])

    def to_dict(self):

        """
        Add stuff for results dictionary here
        """
        result = {'tree': self.tree,
                  'likelihood': self.likelihood}

        return result


class GuenomuParser(object):
    def __init__(self):
        self.best_topology = None
        self.species_trees = None

    def parse(self, filename):
        """
        Extract based on strings defined in __init__
        """
        try:
            self.species_trees = dpy.TreeList.get(path=filename[0],
                                                  schema='nexus')
        except IOError as err:
            logger.error("problem reading species tree file\n" + err)
            return
        try:
            self.best_topology = \
            dpy.TreeList.get(path=filename[1], schema='nexus')[
                0]  # first tree is most frequent (unrooted)
        except IOError as err:
            logger.error("problem reading unrooted species tree file\n" + err)
            return

    def to_dict(self, filenames):
        """
        Add stuff for results dictionary here
        """
        self.parse(filenames)
        result = {
            'best_topology': self.best_topology,
            'species_trees': self.species_trees,
        }
        return result
