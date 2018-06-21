import unittest
import os
import pandas as pd
from rnamake import settings
from rnamake.wrappers import design_rna_wrapper
from rnamake import motif_graph


RESOURCE_PATH = settings.LIB_PATH + "/lib/RNAMake/unittests/unittest_resources/apps/design_rna/"

def _out_files_exist():
    return os.path.isfile("default.out") and os.path.isfile("default.scores")

def _remove_out_files():
    os.remove("default.out")
    os.remove("default.scores")

def _get_motif_graph():
    f = open("default.out")
    lines = f.readlines()
    f.close()

    mg = motif_graph.MotifGraph(lines[0])
    return mg

class DesignRNAUnittest(unittest.TestCase):

    def tests(self):
        path = RESOURCE_PATH + "design_rna_tests.txt"
        df = pd.read_csv(path)

        for i, r in df.iterrows():

            cmd_opts = {
                'pdb'      : RESOURCE_PATH + r.pdb,
                'start_bp' : r.start_bp,
                'end_bp'   : r.end_bp
            }

            print r.pdb

            drw = design_rna_wrapper.DesignRNAWrapper()
            drw.run(**cmd_opts)

            self.failUnless(_out_files_exist(), "out files were not produced")
            mg = _get_motif_graph()

            self.failUnless(len(mg) > 1)
            _remove_out_files()



def main():
    unittest.main()

if __name__ == '__main__':
    main()

