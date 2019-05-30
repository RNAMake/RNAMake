from rnamake import settings
import wrapper

class DDGTectoWrapper(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/ddg_tecto"
        super(self.__class__, self).__init__(program_path)

        self.add_cmd_option("fseq", "", required=True)
        self.add_cmd_option("fss",  "", required=True)
        self.add_cmd_option("cseq", "", required=True)
        self.add_cmd_option("css",  "", required=True)

        self.add_cmd_option("steps", 1000000, required=False)


    def run(self, **options):
        try:
            super(self.__class__, self).run(**options)
        except:
            self._output = ""

    def get_output(self):
        return self._output

    def get_hits(self):
        spl = self._output.split("\n")
        try:
            hits = int(spl[-2])
        except:
            hits = None
        return hits