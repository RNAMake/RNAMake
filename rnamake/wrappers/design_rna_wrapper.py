from rnamake import settings
import wrapper

class DesignRNAWrapper(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/design_rna"
        super(self.__class__, self).__init__(program_path)

        self.add_cmd_option("designs", 1, required=False)
        self.add_cmd_option("seqs_per_design", 1, required=False)
        self.add_cmd_option("out_file", "default.out", required=False)
        self.add_cmd_option("score_file", "default.scores", required=False)
        self.add_cmd_option("verbose", False, required=False)
        self.add_cmd_option("pdbs", False, required=False)

        # no sequence opt
        self.add_cmd_option('only_ideal', False, required=False)

        # start from pdb
        self.add_cmd_option("pdb", "", required=False)
        self.add_cmd_option("start_bp", "", required=False)
        self.add_cmd_option("end_bp", "", required=False)

        self.add_cmd_option("mg", "", required=False)

        # search options
        self.add_cmd_option("search.max_size", 9999, required=False)
        self.add_cmd_option("search.accept_score", 10.0, required=False)


    def get_output(self):
        return self._output



