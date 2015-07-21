import unittest
import random
import rnamake.motif_state_ensemble_tree as motif_state_ensemble_tree
import rnamake.resource_manager as rm
import rnamake.motif_type as motif_type
import rnamake.util as util
import rnamake.basic_io as basic_io
import rnamake.settings as settings
import rnamake.motif as motif
import rnamake.motif_tree as motif_tree
import rnamake.motif_type as motif_type
import rnamake.motif_tree_topology as motif_tree_topology
import build

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
        mst.write_pdbs()

    def test_setup_from_mt(self):
        builder = build.BuildSecondaryStructure()
        ss = builder.build_helix(20)
        con = ss.motif_topology_from_end()
        mtt = motif_tree_topology.MotifTreeTopology(con)
        mt = motif_tree.motif_tree_from_topology_2(mtt)
        mset =  motif_state_ensemble_tree.MotifStateEnsembleTree(mt)

    def test_setup_from_mt_2(self):
        builder = build.BuildMotifTree()
        mt = builder.build()





def main():
    unittest.main()

if __name__ == '__main__':
    main()
