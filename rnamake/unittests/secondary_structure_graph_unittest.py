import unittest
from  rnamake import secondary_structure_graph
import build

class SecondaryStructureGraphUnittest(unittest.TestCase):

    def test_creation(self):
        ssg = secondary_structure_graph.SecondaryStructureGraph()

    def test_from_pose(self):
        builder = build.BuildSecondaryStructure()
        ss_p = builder.build_helix(5)
        ssg = secondary_structure_graph.graph_from_pose(ss_p)





def main():
    unittest.main()

if __name__ == '__main__':
    main()