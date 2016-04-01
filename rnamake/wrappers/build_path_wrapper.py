from rnamake import settings
import wrapper

class BuildPathWrapper(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build_clang/path_builder"
        super(self.__class__, self).__init__(program_path)

        self.add_cmd_option("mg", "", required=True)
