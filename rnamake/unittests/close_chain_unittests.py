import unittest
import numpy as np

from rnamake import chain_closure, resource_manager, motif_graph, atom, util



class ChainClosureUnittests(unittest.TestCase):
    def setUp(self):
        self.rm = resource_manager.ResourceManager()


    def test_simple(self):
        m1 = self.rm.get_bp_step("GG_LL_CC_RR")
        m2 = self.rm.get_bp_step("GU_LL_AC_RR")
        mg = motif_graph.MotifGraph(self.rm)
        mg.add_motif(m1)
        mg.add_motif(m2)
        rna_struct = mg.get_structure()
        #rna_struct.to_pdb("org.pdb")
        p1 = atom.Atom.copy(rna_struct.get_residue(num=508).get_atom("P"))

        chain_closure.close_chain(rna_struct.get_chain(0))
        chain_closure.close_chain(rna_struct.get_chain(1))
        p2 = rna_struct.get_residue(num=508).get_atom("P")
        diff = util.distance(p1.coords, p2.coords)
        self.failUnless(diff > 0.01)
        #rna_struct.to_pdb("test.pdb")


def main():
    unittest.main()

if __name__ == '__main__':
    main()







