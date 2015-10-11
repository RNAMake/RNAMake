import unittest
import build
import rnamake.pose
import rnamake.settings as settings
import rnamake.motif_type as motif_type
import rnamake.motif_tree as motif_tree
import random
import rnamake.eternabot.sequence_designer as sequence_designer
import rnamake.pose_factory as pf


class PoseUnittest(unittest.TestCase):

    def test_creation(self):
        p = pf.factory.pose_from_file("resources/motifs/p4p6")

    def test_designable_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        p = mt.to_pose()
        ss = p.designable_secondary_structure()
        designer = sequence_designer.SequenceDesigner()
        results = designer.design(ss.dot_bracket(), ss.sequence())
        #print results[0]['end'][0]

    def test_motifs(self):
        p = pf.factory.pose_from_file("resources/motifs/p4p6")
        twoways = p.motifs(motif_type.TWOWAY)
        if len(twoways) != 6:
            self.fail("did not properly get all two way junctions")




def main():
    unittest.main()

if __name__ == '__main__':
    main()
