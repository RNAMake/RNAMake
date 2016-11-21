from rnamake import settings
import wrapper

class SimulateTectosWrapper(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/simulate_tectos_devel"
        super(self.__class__, self).__init__(program_path)

        self.add_cmd_option("fseq", "", required=False)
        self.add_cmd_option("fss",  "", required=False)
        self.add_cmd_option("cseq", "", required=False)
        self.add_cmd_option("css",  "", required=False)
        self.add_cmd_option("s", 1000000, required=False)

        self.add_cmd_option("extra_me", "", required=False)

        self.add_cmd_option("new_ggaa_model", False, required=False)
        self.add_cmd_option("ggaa_model", "", required=False)

        # simulation options using for testing
        self.add_cmd_option("simulation.temperature", 298.15, required=False)
        self.add_cmd_option("simulation.steric_radius", 2.2, required=False)
        self.add_cmd_option("simulation.cutoff", 4.5, required=False)


    def run(self, **options):
        try:
            super(self.__class__, self).run(**options)
        except:
            self._output = ""

    def get_output(self):
        spl = self._output.split("\n")
        try:
            hits = int(spl[-2])
        except:
            hits = None
        return hits
