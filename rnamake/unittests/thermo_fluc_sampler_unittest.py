import unittest
import random
import build
import rnamake.thermo_fluc_sampler
import rnamake.motif_tree_topology as motif_tree_topology
import rnamake.motif_state_ensemble_tree as motif_state_ensemble_tree
import rnamake.motif_tree as motif_tree
import rnamake.segmenter as segmenter
import rnamake.resource_manager as rm
import rnamake.pose_factory as pf

class ThermoFlucSamplerUnittest(unittest.TestCase):

    def _fill_basepairs_in_ss(self, ss):
        pairs = ["AU", "UA", "GC", "CG"]
        for bp in ss.basepairs:
            if bp.res1.name == "N" and bp.res2.name == "N":
                p = random.choice(pairs)
                bp.res1.name = p[0]
                bp.res2.name = p[1]

    def test_creation(self):
        tfs = rnamake.thermo_fluc_sampler.ThermoFlucSampler()

    def test_sample(self):
        builder = build.BuildSecondaryStructure()
        ss = builder.build_helix(5)
        con = ss.motif_topology_from_end()
        mtt = motif_tree_topology.MotifTreeTopology(con)
        mt = motif_tree.motif_tree_from_topology(mtt)
        mset =  motif_state_ensemble_tree.MotifStateEnsembleTree(mt)

        tfs = rnamake.thermo_fluc_sampler.ThermoFlucSampler()
        tfs.setup(mset)
        tfs.sample(100)

    def test_relaxer(self):
        s = rnamake.segmenter.Segmenter()
        path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
        p = pf.factory.pose_from_file(path)
        end1 = p.get_basepair(name='A111-A209')[0]
        end2 = p.get_basepair(name='A118-A203')[0]
        segments = s.apply(p, [end1, end2])
        pd = segments.remaining
        pd.name = "p4p6_frag"
        rm.manager.register_motif(pd)
        mt = motif_tree.MotifTree(sterics=0)
        mt.add_motif(pd)
        mt.add_motif(rm.manager.get_motif(name='HELIX.IDEAL.6'),
                     parent_end_name='A111-A209')
        mt.add_connection(0, 1, 'A118-A203')
        ss = mt.designable_secondary_structure()
        self._fill_basepairs_in_ss(ss)
        con = ss.motif_topology_from_end()
        mtt = motif_tree_topology.MotifTreeTopology(con)
        mt2 = motif_tree.motif_tree_from_topology(mtt, sterics=0)

        ni_1 = 0
        ei_1 = 0
        ni_2 = mt2.last_node().index
        ei_2 = 1
        mset =  motif_state_ensemble_tree.MotifStateEnsembleTree(mt2)

        relaxer = rnamake.thermo_fluc_sampler.ThermoFlucRelax()
        relaxer.run(mset, ni_1, ni_2, ei_1, ei_2)
        relaxer.write_pdbs()




def main():
    unittest.main()

if __name__ == '__main__':
    main()