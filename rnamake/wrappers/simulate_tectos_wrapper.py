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
        self.add_cmd_option("temperature", 298.15, required=False)
        self.add_cmd_option("cutoff", 4.5, required=False)
        

        self.add_cmd_option("extra_me", "", required=False)
        self.add_cmd_option("extra_motifs", "", required=False)

        self.add_cmd_option("new_ggaa_model", False, required=False)
        self.add_cmd_option("ggaa_model", "", required=False)

        self.add_cmd_option("record", False, required=False)
        self.add_cmd_option("record_file_type", "", required=False)
        self.add_cmd_option("record_constraints", "", required=False)
        self.add_cmd_option("record_only_bound", False, required=False)
        self.add_cmd_option("record_only_unbound", False, required=False)
        self.add_cmd_option("dump_state", False, required=False)

        self.add_cmd_option("scorer", "", required=False)
        self.add_cmd_option("constraints", "", required=False)



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