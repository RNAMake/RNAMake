import unittest
import os
from rnamake import motif_factory, settings, motif_type


class MotifFactoryUnittest(unittest.TestCase):
    def setUp(self):
        self.mf = motif_factory.MotifFactory()

    def test_creation(self):
        pass

    def test_motifs_from_file(self):
        path = settings.UNITTEST_PATH + "resources/motifs/NWAY.1GID.0"
        motifs = self.mf.motifs_from_file(path)
        self.failUnless(len(motifs) == 3)

        #for i, m in enumerate(motifs):
        #    m.to_pdb("test."+str(i)+".pdb")

    def test_motifs_from_res(self):
        path = settings.UNITTEST_PATH + "resources/motifs/HELIX.IDEAL.6"
        motifs = self.mf.motifs_from_file(path)
        m = motifs[0]

        bps = [ m.get_basepair(0), m.get_basepair(1)]
        res = []

        for bp in bps:
            res.extend(m.get_bp_res(bp))

        new_ms = self.mf.motifs_from_res(res, bps, m, "motif_sub", motif_type.HELIX)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
