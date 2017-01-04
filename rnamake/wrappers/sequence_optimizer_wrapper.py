from rnamake import settings
import wrapper

class SequenceOptimizerWrapper(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/sequence_optimizer_app"
        super(self.__class__, self).__init__(program_path)

        self.add_cmd_option("mg", "", required=True)
        self.add_cmd_option("end_1", "", required=True)
        self.add_cmd_option("end_2", "", required=True)

        self.add_cmd_option("v", False, required=False)
        self.add_cmd_option("out_file", "default.out", required=False)
        self.add_cmd_option("score_file", "default.scores", required=False)
        self.add_cmd_option("n", 1, required=False)
        self.add_cmd_option("pdbs", False, required=False)
        self.add_cmd_option('opt', "Internal", required=False)

        # optimizer options
        self.add_cmd_option("optimizer.cutoff", 10.0, required=False)
        self.add_cmd_option("optimizer.solutions", 1, required=False)
        self.add_cmd_option("optimizer.eterna_cutoff", -1, required=False)
        self.add_cmd_option("optimizer.return_lowest",False, required=False)
        self.add_cmd_option("optimizer.steps", 1000, required=False)


    def run(self, **options):
        try:
            super(self.__class__, self).run(**options)
        except:
            self._output = ""


