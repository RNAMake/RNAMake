import unittest
from rnamake.wrappers import wrapper, build_path_wrapper

class WrapperUnittest(unittest.TestCase):

    def test_creation(self):
        path = "/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/path_builder"
        w = wrapper.Wrapper(path)

    def test_add_cmd_option(self):
        path = "/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/path_builder"
        mg_path = "~/projects/RNAMake.projects/tecto_rna_22_bp/base_mg.top"
        w = wrapper.Wrapper(path)
        w.add_cmd_option("mg", mg_path)
        #print w.get_command()


class BuildPathWrapper(unittest.TestCase):
    def test_creation(self):
        w = build_path_wrapper.BuildPathWrapper()

    def test_get_command(self):
        w = build_path_wrapper.BuildPathWrapper()

        try:
            w.get_command()
            raise RuntimeError
        except wrapper.WrapperException:
            pass
        except:
            self.fail("did not get correct error")

        w.set_cmd_option("mg", "~/projects/RNAMake.projects/tecto_rna_22_bp/base_mg.top")
        cmd = w.get_command()
        if cmd != "./Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/path_builder -mg \"~/projects/RNAMake.projects/tecto_rna_22_bp/base_mg.top\"":
            print "made it"

    def _test_run(self):
        w = build_path_wrapper.BuildPathWrapper()
        mg_path = "/Users/josephyesselman/projects/RNAMake.projects/tecto_rna_22_bp/base_mg.top"
        w.run(mg=mg_path)




def main():
    unittest.main()

if __name__ == '__main__':
    main()