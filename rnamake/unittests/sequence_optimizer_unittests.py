import unittest

from rnamake import motif_graph, motif_topology
from rnamake import sequence_optimizer
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

    def test_init_2(self):
        builder = build.BuildMotifGraph()
        mg = builder.build(3)
        mg.replace_ideal_helices()

        c = motif_topology.GraphtoTree()
        mt = c.convert(mg)
        build_points = mt.get_build_points()

        so = sequence_optimizer.SequenceOptimizer3D()
        so.get_optimized_sequences(mt, build_points[0].node.data.ends[1])


def main():
    unittest.main()

if __name__ == '__main__':
    main()
