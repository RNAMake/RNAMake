import unittest

from rnamake import motif_graph, motif_topology
from rnamake import sequence_optimizer, settings
from rnamake import resource_manager as rm
import build

class SequenceOptimizerUnittests(unittest.TestCase):

    def setUp(self):
        pass

    def _test_init(self):
        mg = motif_graph.MotifGraph()
        mg.add_motif(m_name="HELIX.IDEAL.10")
        mg.replace_ideal_helices()

        c = motif_topology.GraphtoTree()
        mt = c.convert(mg)
        build_points = mt.get_build_points()

        so = sequence_optimizer.SequenceOptimizer3D()
        so.get_optimized_sequences(mt, build_points[0].node.data.ends[1])

    def _test_init_2(self):
        builder = build.BuildMotifGraph()
        mg = builder.build(5)
        mg.add_motif(m_name="HAIRPIN.1C0A.0")
        mg.replace_ideal_helices()
        c = motif_topology.GraphtoTree()
        mt = c.convert(mg)
        mt.write_pdbs()

        so = sequence_optimizer.SequenceOptimizer3D()
        solutions = so.get_optimized_sequences(mt,mt.last_node().data.ends[0],
                                               mt.last_node().index, 0)

    def _test_minittr(self):
        path = settings.UNITTEST_PATH + "/test_problems/mini_ttr/"
        f = open(path+"sol.mg")
        lines = f.readlines()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=lines[0])
        mg.replace_ideal_helices()
        n1 = mg.get_node(1)
        n2 = n1.connections[1].partner(n1.index)

        c = motif_topology.GraphtoTree()
        mt = c.convert(mg, last_node=n2)
        mt.option('sterics', 0)

        so = sequence_optimizer.SequenceOptimizer3D()


        solutions = so.get_optimized_sequences(mt,mt.last_node().data.ends[0],
                                               mt.get_node(4).index, 1)

        print len(solutions)

    def test_minittr_2(self):
        path = "/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build"
        f = open(path+"/default.out.bak")
        lines = f.readlines()
        f.close()

        mg = motif_graph.MotifGraph(mg_str=lines[0])

        n1 = mg.get_node(0)
        n2 = n1.connections[0].partner(n1.index)

        c = motif_topology.GraphtoTree()
        mt = c.convert(mg, last_node=n2)
        mt.option('sterics', 0)

        so = sequence_optimizer.SequenceOptimizer3D()

        solutions = so.get_optimized_sequences(mt,mt.get_node(0).data.ends[0],
                                               mt.last_node().index, 1)

        dss = mg.designable_secondary_structure()
        dss.replace_sequence(solutions[0].sequence)
        mg.replace_helix_sequence(dss)
        mg.write_pdbs()


def main():
    unittest.main()

if __name__ == '__main__':
    main()
