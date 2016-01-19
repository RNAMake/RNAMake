import unittest
import build
from rnamake import secondary_structure_tree, motif_tree

class SecondaryStructureTreeUnittest(unittest.TestCase):

    def test_creation(self):
        sst = secondary_structure_tree.SecondaryStructureTree()

    def test_from_pose(self):
        builder = build.BuildSecondaryStructure()
        ss_p = builder.build_helix(10)
        sst = secondary_structure_tree.tree_from_pose(ss_p)
        if len(ss_p.motifs) != len(sst.tree):
            self.fail("did not get the right number of motifs in tree")

    def test_convert_to_motifs(self):
        builder = build.BuildSecondaryStructure()
        ss_p = builder.build_helix(10)
        sst = secondary_structure_tree.tree_from_pose(ss_p)
        mt = motif_tree.motif_tree_from_ss_tree(sst)
        if len(mt) != 9:
            self.fail("did not get the correct number of nodes in motif tree")

def main():
    unittest.main()

if __name__ == '__main__':
    main()