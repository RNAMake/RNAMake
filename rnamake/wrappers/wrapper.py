import os
import subprocess
from options import Options, OptionType

class WrapperException(Exception):
    pass

class Wrapper(object):

    __slots__ = ['__name',
                 '__program_path',
                 '__options',
                 '__cmd_options',
                 '__default_cmd_options',
                 '__required_cmd_options']

    def __init__(self, program_path, **cmd_options):
        self.__name = "Wrapper"
        self.__program_path = program_path
        if not os.path.exists(self.__program_path):
            raise ValueError("program path : " + program_path + " does not exist")
        self.__cmd_options = Options()
        self.__default_cmd_options = {}
        self.__required_cmd_options = {}

        self.__options = Options()
        self.__options.add("ignore_defaults", True)

    def add_cmd_option(self, name, value, required=False):
        self.__cmd_options.add(name, value)
        self.__default_cmd_options[name] = value
        self.__required_cmd_options[name] = required

    def set_cmd_option(self, name, value):
        self.__cmd_options[name] = value

    def get_command(self, **options):
        self.__cmd_options.dict_set(options)

        s = self.__program_path + " "
        for k, opt in self.__cmd_options:
            if self.__is_required_option(k) and self.__is_default_option(k):
                if self.__default_cmd_options[k] == opt.value:
                    raise WrapperException("required option: " + k + " was not supplied")

            s += "-" + k + " "
            if opt.otype == OptionType.STRING:
                s += "\"" + opt.value + "\" "
            if opt.otype == OptionType.NUMBER:
                s += str(opt.value) + " "
        return s

    def run(self, **options):
        cmd = self.get_command(**options)
        subprocess.call(cmd, shell=True)

    def __is_required_option(self, name):
        return self.__required_cmd_options[name]

    def __is_default_option(self, name):
        if name in self.__default_cmd_options:
            return 1
        else:
            return 0


