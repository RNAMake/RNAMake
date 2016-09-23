import unittest
import os
import random
import rnamake.resource_manager as rm
import rnamake.motif_tree as motif_tree

from rnamake import sqlite_library
from rnamake import motif_graph, motif_topology, motif_state_ensemble_tree

class MotifStateEnsembleTreeUnittest(unittest.TestCase):

    def test_creation(self):
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree()

    def test_add(self):
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree()
        mse = rm.manager.get_motif_state_ensemble(name="GG_LL_CC_RR")
        mset.add_ensemble(mse)
        mset.add_ensemble(mse)
        if len(mset) != 2:
            self.fail("did not properly add ensembles")

    def test_to_mst(self):
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree()
        mse = rm.manager.get_motif_state_ensemble(name="GG_LL_CC_RR")
        mset.add_ensemble(mse)
        mset.add_ensemble(mse)
        mset.add_ensemble(mse)
        mset.add_ensemble(mse)

        mst = mset.to_mst()
        #print mst.get_node(0).data.name()
        #mst.write_pdbs()

    def _test_enumerator(self):
        lib = sqlite_library.MotifStateEnsembleSqliteLibrary("all_bp_steps")
        mtst = motif_state_ensemble_tree.MotifStateEnsembleTree()
        mtst.add_ensemble(lib.get(name="GG_LL_CC_RR"))
        mtst.add_ensemble(lib.get(name="GG_LL_CC_RR"))
        enumerator = motif_state_ensemble_tree.MotifStateEnsembleTreeEnumerator(mtst)
        enumerator.record()

    def test_setup_from_mt_4(self):
        mg = motif_graph.MotifGraph()
        mg.add_motif(m_name="HELIX.IDEAL.20")
        mg.replace_ideal_helices()
        mt = motif_topology.graph_to_tree(mg)
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt)

        mst = mset.to_mst()
        #for r in mt.residues():
        #    print mst.get_residue(r.uuid)

    def test_motif_sub(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        m1 = mlib.get_random()
        m2 = mlib.get_random()
        while len(m1.sequence()) == len(m2.sequence()):
            m2 = mlib.get_random()

        motif_state_ensemble_tree.build_motif_sub_for_motif_state_ensemble(m1, [m2], mlib)
        rm.manager.register_extra_motif_ensembles("test.dat")

        r = rm.manager.has_supplied_motif_ensemble(m1.name, m1.ends[0].name())
        self.failUnless(r == 1, "did not find sub")

        r = rm.manager.has_supplied_motif_ensemble(m1.name, m1.ends[1].name())
        self.failUnless(r == 1, "did not find sub")

        mt = motif_tree.MotifTree()
        mt.add_motif(m_name="HELIX.IDEAL")
        mt.add_motif(m1)
        mt.add_motif(m_name="HELIX.IDEAL")

        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt=mt)
        mst = mset.to_mst()

        self.failUnless(mst.get_node(1).data.cur_state.name == m2.name)

        os.remove("test.dat")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
