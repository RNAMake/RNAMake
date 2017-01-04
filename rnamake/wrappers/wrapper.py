import os
import subprocess
from options import Options, OptionType

class WrapperException(Exception):
    pass

class Wrapper(object):

    __slots__ = ['_name',
                 '_program_path',
                 '_options',
                 '_output',
                 '_cmd_options',
                 '_default_cmd_options',
                 '_required_cmd_options']

    def __init__(self, program_path, **cmd_options):
        self._name = "Wrapper"
        self._output = ""
        self._program_path = program_path
        if not os.path.exists(self._program_path):
            raise ValueError("program path : " + program_path + " does not exist")
        self._cmd_options = Options()
        self._default_cmd_options = {}
        self._required_cmd_options = {}

        self._options = Options()
        self._options.add("ignore_defaults", True)

    def setup_from_file(self, fname):
        f = open(fname)
        lines = f.readlines()
        f.close()
        for l in lines:
            spl = l.split()
            cmd_name = spl[0][1:]
            if len(spl) == 2:
                val = spl[1]
                org_val = self.get_cmd_option(cmd_name)
                if type(org_val) == int:
                    val = int(spl[1])
                elif type(org_val) == float:
                    val = float(spl[1])
                self.set_cmd_option(cmd_name, val)
            else:
                self.set_cmd_option(cmd_name, True)

    def add_cmd_option(self, name, value, required=False):
        self._cmd_options.add(name, value)
        self._default_cmd_options[name] = value
        self._required_cmd_options[name] = required

    def get_cmd_option(self, name):
        return self._cmd_options[name]

    def set_cmd_option(self, name, value):
        self._cmd_options[name] = value

    def get_command(self, **options):
        self._cmd_options.dict_set(options)

        s = self._program_path + " "
        for k, opt in self._cmd_options:
            if self._is_required_option(k) and self._is_default_option(k):
                if self._default_cmd_options[k] == opt.value:
                    raise WrapperException("required option: " + k + " was not supplied")

            if self._options['ignore_defaults'] and self._is_default_option(k) and \
               opt.value == self._default_cmd_options[k]:
                continue

            s += "-" + k + " "
            if opt.otype == OptionType.STRING:
                s += "\"" + opt.value + "\" "
            if opt.otype == OptionType.NUMBER:
                s += str(opt.value) + " "
        return s

    def run(self, **options):
        cmd = self.get_command(**options)
        self._output = subprocess.check_output(cmd, shell=True)
        return self._output

    def _is_required_option(self, name):
        return self._required_cmd_options[name]

    def _is_default_option(self, name):
        if name in self._default_cmd_options:
            return 1
        else:
            return 0


