from rnamake import settings
import wrapper

class SimulateTectosWrapper(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/simulate_tectos_devel"
        super(self.__class__, self).__init__(program_path)

