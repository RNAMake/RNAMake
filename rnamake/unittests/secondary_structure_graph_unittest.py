import unittest
from  rnamake import secondary_structure_graph
import build

class SecondaryStructureGraphUnittest(unittest.TestCase):

    def test_creation(self):
        ssg = secondary_structure_graph.SecondaryStructureGraph()

    def test_from_pose(self):
        builder = build.BuildSecondaryStructurePose()
        ss_p = builder.build_helix(5)
        print len(ss_p.motifs)
        ssg = secondary_structure_graph.graph_from_pose(ss_p)
        print ssg





def main():
    unittest.main()

if __name__ == '__main__':
    main()