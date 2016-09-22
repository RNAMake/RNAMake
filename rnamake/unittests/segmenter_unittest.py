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
        path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
        p = pf.factory.pose_from_file(path)
        twoways = p.motifs(motif_type.TWOWAY)
        for i, t in enumerate(twoways):
            s = rnamake.segmenter.Segmenter()
            segments = s.apply(p, t.ends)
            self.failUnless(len(t.residues()) == len(segments.removed.residues()))


def main():
    unittest.main()

if __name__ == '__main__':
    main()
