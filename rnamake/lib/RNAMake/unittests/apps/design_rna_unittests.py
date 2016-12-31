import unittest
import os
from rnamake import settings
from rnamake.wrappers import design_rna_wrapper
from rnamake import motif_graph


#./design_rna -pdb ~/Downloads/3suy.pdb -end_bp "X30-X56" -start_bp "X20-X78" -pdbs -verbose -search.max_size 90

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

    def test_1(self):

        cmd_opts = {
            'pdb'      : RESOURCE_PATH + '3suy.pdb',
            'start_bp' : 'X20-X78',
            'end_bp'   : 'X30-X56',
            'search.max_size' : 90}

        drw = design_rna_wrapper.DesignRNAWrapper()
        drw.run(**cmd_opts)

        self.failUnless(_out_files_exist(), "out files were not produced")
        mg = _get_motif_graph()

        self.failUnless(len(mg) > 1)

        _remove_out_files()

    def test_2(self):
        cmd_opts = {
            'pdb'      : RESOURCE_PATH + '1xpe.pdb',
            'start_bp' : 'B5-B19',
            'end_bp'   : 'A5-A19',
            'search.max_size' : 50}

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

