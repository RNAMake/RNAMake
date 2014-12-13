import abc

class Base(object):
    def __init__(self, **options):
        pass

    @abc.abstractmethod
    def setup_options_and_constraints(self):
        self.options = self.set_default_options()
        self.constraints = {}

    @abc.abstractmethod
    def set_default_options(self):
        return { 'test' : None }

    def set_options(self, options, error=1):
        for key in options:
            if   key in self.options:
                self.options[key] = options[key]
            elif not error:
                pass
            else:
                raise ValueError("invalid option " + key)

    def option(self, option, value=None):
        if key in self.options:
            if value is None:
                return self.options[option]
            else:
                self.options[option] = key
        else:
            raise ValueError("invalid option " + option)






