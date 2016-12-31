import unittest

from rnamake import motif_graph, motif_topology, motif_state_graph
from rnamake import sequence_optimizer, settings, util
import is_equal
from rnamake import resource_manager as rm
import build

class SequenceOptimizerUnittests(unittest.TestCase):

    def setUp(self):
        pass

    def test_simple_hairpin(self):
        mg = motif_graph.MotifGraph()
        for i in range(10):
            mg.add_motif(m_name="HELIX.IDEAL")
        mg.add_motif(m_name="HAIRPIN.1GID.0")

        so = sequence_optimizer.SequenceOptimizer3D()
        scorer = sequence_optimizer.ExternalTargetScorer(
                    mg.get_node(9).data.ends[1].state().copy(), 8, 1)
        sols = so.get_optimized_sequences(mg, scorer)
        bp1 = mg.get_node(9).data.ends[1].copy()
        mg.replace_helix_sequence(seq=sols[0].sequence)

        self.failUnless(mg.sequence() == sols[0].sequence)
        bp2 = mg.get_node(8).data.ends[1]
        self.failUnless(abs(sols[0].dist_score - bp1.diff(bp2)) < 0.2)

    def test_simple_hairpin_2(self):
        mg = motif_graph.MotifGraph()
        for i in range(10):
            mg.add_motif(m_name="HELIX.IDEAL")
        mg.add_motif(m_name="HAIRPIN.1GID.0")
        so = sequence_optimizer.SequenceOptimizer3D()
        scorer = sequence_optimizer.ExternalTargetScorer(
                    mg.get_node(9).data.ends[1].state().copy(), 8, 1)
        mg_opt = so.get_optimized_mg(mg, scorer)

        mg.replace_helix_sequence(seq=mg_opt.sequence())

        for n in mg:
            self.failUnless(n.data.name == mg_opt.get_node(n.index).data.name)

        atoms1 = mg.get_structure().structure.atoms()
        atoms2 = mg_opt.get_structure().structure.atoms()

        for i in range(len(atoms1)):
            diff = util.distance(atoms1[i].coords, atoms2[i].coords)
            self.failUnless(diff < 0.1)

    def test_minittr(self):
        path = settings.UNITTEST_PATH + "/test_problems/mini_ttr/"
        f = open(path+"sol.mg")
        lines = f.readlines()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=lines[0])
        mg.replace_ideal_helices()
        scorer = sequence_optimizer.InternalTargetScorer(11, 1, 19, 1)
        so = sequence_optimizer.SequenceOptimizer3D()
        mg_opt = so.get_optimized_mg(mg, scorer)

        mg.replace_helix_sequence(seq=mg_opt.sequence())

        for n in mg:
            self.failUnless(n.data.name == mg_opt.get_node(n.index).data.name)

        atoms1 = mg.get_structure().structure.atoms()
        atoms2 = mg_opt.get_structure().structure.atoms()

        for i in range(len(atoms1)):
            diff = util.distance(atoms1[i].coords, atoms2[i].coords)
            self.failUnless(diff < 0.1)

    def test_chip_only(self):
        base_dir = settings.UNITTEST_PATH + "resources/motif_graph/"
        f = open(base_dir+"tecto_chip_only.mg")
        l = f.readline()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=l)
        #for n in mg:
        #    print n.index, n.data.name

        scorer = sequence_optimizer.InternalTargetScorer(22, 1, 19, 1)
        so = sequence_optimizer.SequenceOptimizer3D()
        so.option('return_lowest', 1)
        so.option('max_steps', 100)
        mg_opt = so.get_optimized_mg(mg, scorer)

        mg.replace_helix_sequence(seq=mg_opt.sequence())

        mg.write_pdbs()
        mg_opt.write_pdbs("opt")

        self.failUnless(mg_opt.sequence() == mg.sequence())

        for n in mg:
            self.failUnless(n.data.name == mg_opt.get_node(n.index).data.name)




def main():
    unittest.main()

if __name__ == '__main__':
    main()
