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
    __slots__ = ["value", "otype"]

    def __init__(self, value, otype):
        self.value, self.otype = value, otype

class Options(object):
    """
    holds the options for a class derived from Base. Checks for correct type
    of option values

    Attributes
    ----------
    `options` : Dict of Option objects
        stores all the current options
    """
    def __init__(self, options={}):
        self.options = {}
        self.dict_add(options)

    def __repr__(self):
        pass

    def dict_add(self, options):
        for k,v in options.iteritems():
            self.add(k, v)

    def dict_set(self, options):
        for k,v in options.iteritems():
            self.set(k, v)

    def add(self, name, value):
        if name in self.options:
            raise ValueError("cannot add option "+ name +", already exists")
        otype = type(value)

        self.options[name] = Option(value, otype)

    def set(self, name, value):
        if name not in self.options:
            raise ValueError("cannot set option "+ name +", does not exists")
        option = self.options[name]

        if option.otype != type(value) :
            raise ValueError("cannot set option "+ str(name) +", it is typed " +
                             "as an " + str(option.otype) + ", but was " +
                             "supplied " + str(type(value)))

        self.options[name].value = value

    def get(self, name):
        if name not in self.options:
            raise ValueError("cannot get option "+ name +", it does not exist")

        return self.options[name].value

    def valid_options(self):
        return self.options.keys()

    def __contains__(self, name):
        return name in self.options
