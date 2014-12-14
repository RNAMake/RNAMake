import abc
from . import option

class Base(object):
    def __init__(self, **options):
        self.setup_options_and_constraints()

    @abc.abstractmethod
    def setup_options_and_constraints(self):
        self.options = self.set_default_options()
        self.constraints = {}

    @abc.abstractmethod
    def set_default_options(self):
        return { 'test' : None, 'test2' : None }

    def set_options(self, options, error=1):
        for key in options:
            if   key in self.options:
                self.options[key] = options[keyi]
            elif not error:
                pass
            else:
                raise ValueError("invalid option " + key)

    def set_constraints(self, constraints, error=1):
        for key in constraints:
            if   key in self.constraints:
                self.constraints[key] = constraints[key]
            elif not error:
                pass
            else:
                raise ValueError("invalid constraint " + key)

    def option(self, option, value=None):
        if option in self.options:
            if value is None:
                return self.options[option]
            else:
                self.options[option] = value
        else:
            raise ValueError("invalid option " + option)






