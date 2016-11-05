class OptionType(object):
    NUMBER = 0
    STRING = 1
    BOOL   = 2

class OptionException(Exception):
    pass


class Option(object):
    """
    stores the an option for a class dervived from Base. contains both its
    value and type(value) to ensure that an option is not set to some
    bizzare value. Option objects should not be generated directly but
    will be made in an Options object that each class dervived from Base
    should contain

    :param value: the value of a given option to be stored
    :param otype: the type of the value

    :type value: unknown
    :type otype: Type object

    Attributes
    ----------
    `value` : Unknown type
        the value of the option to be stored
    `otype` : Type object
        the type of a given value
    """

    __slots__ = ["__name", "__value", "__otype"]

    def __init__(self, value, otype):
        self.__otype =  otype
        # explicitily calling @value.setter to make sure value is of correct
        # type
        self.value = value

    @property
    def value(self):
        return self.__value

    @property
    def otype(self):
        return self.__otype

    @value.setter
    def value(self, new_value):

        if self.__otype == OptionType.NUMBER and \
           not (isinstance(new_value, int) or isinstance(new_value, float)):
            raise OptionException("attemped to set value to incorrect type")

        #True and False are counted as ints for some reason
        #if self.__otype == OptionType.NUMBER and \
        #   (new_value == True or new_value == False):
        #    raise OptionException("attemped to set value to incorrect type")

        if self.__otype == OptionType.BOOL and \
           not isinstance(new_value, bool):
            raise OptionException("attemped to set value to incorrect type")

        if self.__otype == OptionType.STRING and \
           not isinstance(new_value, str):
            raise OptionException("attemped to set value to incorrect type")

        self.__value = new_value



class Options(object):
    """
    holds the options for a class derived from Base. Checks for correct type
    of option values

    Attributes
    ----------
    `options` : Dict of Option objects
        stores all the current options
    """

    __slots__ = ["__options"]

    def __init__(self, options={}):
        self.__options = {}
        self.dict_add(options)

    def __repr__(self):
        pass

    def dict_add(self, options):
        for k,v in options.iteritems():
            self.add(k, v)

    def dict_set(self, options, error=1):
        for k,v in options.iteritems():
            try:
                self[k] = v
            except OptionException as e :
                if error:
                    raise e

    def add(self, name, value):
        if name in self.__options:
            raise OptionException("cannot add option "+ name +", already exists")

        otype = None
        if isinstance(value, int) or isinstance(value, float):
            otype = OptionType.NUMBER
        if isinstance(value, str):
            otype = OptionType.STRING
        if isinstance(value, bool):
            otype = OptionType.BOOL

        if otype is None:
            raise OptionException("not a valid option value")

        self.__options[name] = Option(value, otype)

    def valid_options(self):
        return self.__options.keys()

    def get_dict(self):
        opts = {}
        for name, opt in self.__options.iteritems():
            opts[name] = opt.value
        return opts

    def __contains__(self, name):
        return name in self.__options

    def __getitem__(self, name):
        if name not in self.__options:
            raise OptionException("cannot get option "+ name + ", it does not exist")

        return self.__options[name].value

    def __setitem__(self, name, value):
        if name not in self:
            raise OptionException("cannot set option "+ name +", does not exists")

        self.__options[name].value = value

    def __iter__(self):
        return self.__options.iteritems()