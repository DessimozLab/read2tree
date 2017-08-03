from numbers import Integral, Real
from six import string_types
from abc import ABCMeta, abstractproperty
from dendropy import Tree


class Option(object):
    """Abstract base class for an option.

    Options provide an interface between the wrapper and the
    concrete command line option of the wrapped program."""
    __metaclass__ = ABCMeta

    def __init__(self, name, default=None, active=False):
        self._name = name
        self.set_value(default)
        self.active = active

    def __repr__(self):
        return '{}({}={}) <{}>'.format(self.__class__.__name__, self.name, self.get_value(), 'on' if self.active else 'off')

    def __str__(self):
        return (' '.join([self._name, str(self.get_value())]) if self.active else '')

    @property
    def active(self):
        return self._active

    @active.setter
    def active(self, val):
        self._active = True if val else False

    @property
    def name(self):
        return self._name

    def set_value(self, value):
        self._value = value
        if value is not None:
            self.active = True

    def get_value(self):
        return self._value

    def set_and_activate(self, value):
        self.set_value(value)
        self.active = True

    def status(self):
        return 'Name: {}\nValue: {}\nActive: {}\nStr: {}'.format(self.name,
                                                                 self.get_value(),
                                                                 self.active,
                                                                 str(self) or "''")


class ValueOption(Option):
    __metaclass__ = ABCMeta


class TypedValueOption(ValueOption):
    """A TypedValueOption is an option that only accepts options of a given type.

    This abstract class provides the functionality to check the type
    of a passed value and raises an ValueError if it doesn't match
    the expected type.

    A TypedValueOption must overwrite the abstract property _type.
    """

    __metaclass__ = ABCMeta

    @abstractproperty
    def _type(self):
        pass

    def set_value(self, value):
        if isinstance(value, self._type):
            self._value = value
            self.active = True

        else:
            raise ValueError('Value should be of type {}'.format(self.type))


### Concrete classes from here on

class IntegerOption(TypedValueOption):
    """option to hold an integer value"""
    @property
    def _type(self):
        return Integral


class FloatOption(TypedValueOption):
    """Option to hold a real number value"""

    @property
    def _type(self):
        return Real

    def get_value(self):
        return float(self._value)


class StringOption(TypedValueOption):
    """Opion to hold a string value"""

    def __init__(self, name, value=None, active=False):
        if value is None:
            value = str()
        super(StringOption, self).__init__(name, value, active)

    @property
    def _type(self):
        return string_types


class FlagOption(TypedValueOption):
    """Option to hold a boolean flag value, i.e. True or False"""
    @property
    def _type(self):
        return bool

    def __str__(self):
        return (self._name if self.active and self.get_value() else '')

    def set_value(self, value):
        self._value = bool(value)

    def get_value(self):
        return self._value


class TreeInputOption(TypedValueOption):
    """Option to hold a phylogenetic tree argument.

    As of now, Trees are represented as :class:`dendropy.Tree` objects."""

    @property
    def _type(self):
        return Tree


class MultiOption(Option):
    """Option to hold a list"""

    @property
    def _type(self):
        return list

    def __str__(self):
        listopts = self.get_value()
        if listopts is None: return ''
        strings = []
        for item in listopts:
            item_string = ' '.join([self._name, str(item)]) if self.active else ''
            if item_string > '':
                strings.append(item_string)

        return ' '.join(strings)


class OptionSet(object):
    """Option to hold a set of key-value pairs."""
    def __init__(self, options):
        if isinstance(options, (list, tuple)):
            self.options = {opt.name: opt for opt in options}
        elif isinstance(options, dict):
            self.options = options
        else:
            raise ValueError('Expected a list, tuple or dict of options, not {}'.format(type(options)))

    def __str__(self):
        strings = []
        for name, option in self.options.items():
            option_string = str(option)
            if option_string > '':
                strings.append(option_string)

        return ' '.join(strings)

    def __getitem__(self, item):
        return self.options[item]

    def list(self):
        return [(name, option) for (name, option) in self.options.items()]
