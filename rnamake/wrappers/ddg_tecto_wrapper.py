from rnamake import settings
import wrapper

from collections import namedtuple

class CPPOption(object):
    def __init__(self, name, default_val, type, required):
        self.name = name
        self.default_val = default_val
        self.type = type
        self.required = required

def parse_option_type(type_name):
    if   type_name == "base::OptionType::STRING":
        return wrapper.OptionType.STRING
    elif type_name == "base::OptionType::INT":
        return wrapper.OptionType.NUMBER
    elif type_name == "base::OptionType::FLOAT":
        return wrapper.OptionType.NUMBER
    elif type_name == "base::OptionType::BOOL":
        return wrapper.OptionType.BOOL
    else:
        raise ValueError("unknown type: " + type_name)

def parse_options_from_cpp(path):
    f = open(path)
    lines = f.readlines()
    f.close()

    cpp_options = []

    for l in lines:
        if l.find("add_option") != -1:
            spl = l[15:-3].split(",")
            name = spl[0].rstrip().lstrip()[1:-1]
            type = parse_option_type(spl[2].rstrip().lstrip())
            default_val = ""

            if   type == wrapper.OptionType.STRING:
                default_val = spl[1].rstrip().lstrip()[1:-1]
            elif type == wrapper.OptionType.NUMBER:
                default_val = float(spl[1])
            elif type == wrapper.OptionType.BOOL:
                default_val = bool(spl[2])

            required = False
            if len(spl) == 4:
                if spl[-1].rstrip().lstrip() == 'true':
                    required = True

            cpp_options.append(CPPOption(name, default_val, type, required))

    return cpp_options

DDGTectoResults = namedtuple('DDGTectoResults', ['runs', 'avg', 'stdev'])

class DDGTectoWrapper(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/ddg_tecto"
        super(self.__class__, self).__init__(program_path)
        self._error = ""

        path_to_cpp = "/lib/RNAMake/apps/ddg_tecto/ddg_tecto.cc"
        cpp_options = parse_options_from_cpp(settings.LIB_PATH + path_to_cpp)

        for cpp_opt in cpp_options:
            self.add_cmd_option(cpp_opt.name, cpp_opt.default_val,
                                required=cpp_opt.required)

    def run(self, **options):
        super(self.__class__, self).run(**options)

    def get_output(self):
        return self._output

    def get_errors(self):
        spl = self._output.split("\n")
        for l in spl:
            if l.find('print_bracktrace') != -1:
                self._error = l
                break

        if len(self._error) > 0:
            return self._error

    def get_results(self):
        self.get_errors()
        if len(self._error) > 0:
            raise wrapper.WrapperException(
                "cannot return results an errror was returned: " + self._error)

        spl = self._output.split("\n")
        result_line = spl[-2].split()
        return DDGTectoResults(int(result_line[4]), float(result_line[6]),
                              float(result_line[8]))


def main():
    dtw = DDGTectoWrapper()


if __name__ == "__main__":
    main()





































