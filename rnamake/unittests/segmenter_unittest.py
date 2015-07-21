import unittest
import rnamake.segmenter
import rnamake.pose
import rnamake.settings
import rnamake.pose_factory as pf
import rnamake.motif_type as motif_type

class SegmenterUnittest(unittest.TestCase):

    def test_creation(self):
        s = rnamake.segmenter.Segmenter()

    def test_apply(self):
        s = rnamake.segmenter.Segmenter()
        path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
        p = pf.factory.pose_from_file(path)
        nways = p.motifs(motif_type.NWAY)
        for t in nways:
            segments = s.apply(p, t.ends)
            segments.remaining.to_pdb("remaining.pdb")
            segments.removed.to_pdb("removed.pdb")



def main():
    unittest.main()

if __name__ == '__main__':
    main()
