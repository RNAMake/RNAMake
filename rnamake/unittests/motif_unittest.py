import unittest
import os
import rnamake.motif
import rnamake.settings
import util

class MotifUnittest(unittest.TestCase):

    def setUp(self):
        path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
        self.motif = util.supress_log_output(rnamake.motif.Motif, path)

    def test_creation(self):
        m = rnamake.motif.Motif()

    def test_creation_mdir(self):
        try:
            path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
            m = rnamake.motif.Motif(path)
        except:
            self.fail("did not generate motif correctly")

    def test_create_pdb(self):
        try:
            path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
            m = rnamake.motif.Motif(pdb=path)
            os.remove("p4p6_dssr.out")
            os.remove("ref_frames.dat")
        except:
            self.fail("did not generate motif correctly")

    def test_get_basepair_ends(self):
         m = self.motif
         m.setup_basepair_ends()

    def test_get_basepair(self):
        m = self.motif
        bp = m.basepairs[0]

        found = m.get_basepair(res1=bp.res1, res2=bp.res2)
        print len(found)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
