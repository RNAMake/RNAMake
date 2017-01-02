import unittest
import os
from rnamake import settings
from rnamake.wrappers import wrapper, simulate_tectos_wrapper

class WrapperUnittest(unittest.TestCase):
    def setUp(self):
        path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/design_rna"
        if not os.path.isfile(path):
            self.skipTest("executable not available")

    def test_creation(self):
        path = "/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/design_rna"
        w = wrapper.Wrapper(path)

    def test_add_cmd_option(self):
        path = "/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/design_rna"
        mg_path = "~/projects/RNAMake.projects/tecto_rna_22_bp/base_mg.top"
        w = wrapper.Wrapper(path)
        w.add_cmd_option("mg", mg_path)
        #print w.get_command()

    def test_command(self):
        path = "/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/design_rna"
        w = wrapper.Wrapper(path)
        #cmd = w.get_command(fake_arg=2)
        #print cmd


class BuildPathWrapper(unittest.TestCase):
    def setUp(self):
        path = settings.LIB_PATH + "/rnamake/lib/RNAMake/cmake/build/path_builder"
        if not os.path.isfile(path):
            self.skipTest("executable not available")

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
            pass


class SimulateTectosWrapper(unittest.TestCase):
    def setUp(self):
        path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/simulate_tectos_devel"
        if not os.path.isfile(path):
            self.skipTest("executable not available")

    def test_creation(self):
        w = simulate_tectos_wrapper.SimulateTectosWrapper()
        
    def _test_run(self):
        w = simulate_tectos_wrapper.SimulateTectosWrapper()
        w.run()
        hits = w.get_output()
        self.failUnless(hits > 0)

    def test_with_options(self):
        # really long tecto expect no hits, make sure its getting new sequence
        seq = 'CUAGGAUAUGGAAGUGGGCUUCGGGAACGAAGCCCACUUCCUAAGUCCUAG'
        ss  = '(((((((..((((((((((((((....))))))))))))))...)))))))'
        steps = 1000

        w = simulate_tectos_wrapper.SimulateTectosWrapper()
        w.run(cseq=seq, css=ss, s=steps)
        hits = w.get_output()
        self.failUnless(hits < 5)

        # test with dictionary instead
        opt_dict = {'cseq' : 'CUAGGAUAUGGAAGUGGGCUUCGGGAACGAAGCCCACUUCCUAAGUCCUAG',
                    'css'  : '(((((((..((((((((((((((....))))))))))))))...)))))))',
                    's'    : 1000 }

        w = simulate_tectos_wrapper.SimulateTectosWrapper()
        w.run(**opt_dict)
        hits = w.get_output()
        self.failUnless(hits < 5)



    def test_with_options_from_dict(self):
        try:
            import pandas as pd
        except:
            self.skipTest("pandas is not available ")

        df = pd.DataFrame(columns="cseq css s".split())
        df.loc[0] = ['CUAGGAUAUGGAAGUGGGCUUCGGGAACGAAGCCCACUUCCUAAGUCCUAG',
                     '(((((((..((((((((((((((....))))))))))))))...)))))))',
                     1000]

        w = simulate_tectos_wrapper.SimulateTectosWrapper()
        w.run(**df.loc[0].to_dict())
        hits = w.get_output()
        self.failUnless(hits < 5)


def main():
    unittest.main()

if __name__ == '__main__':
    main()