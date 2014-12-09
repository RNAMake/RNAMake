import unittest
import os
import rnamake.x3dna

class X3dnaUnittest(unittest.TestCase):

    def test_generate_ref_frame(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6"
        x3dna = rnamake.x3dna.X3dna()
        x3dna.generate_ref_frame(path)
        if not os.path.isfile("ref_frames.dat"):
            self.fail("ref_frames.dat file should of been generated")
        if os.path.isfile("bestpairs.pdb"):
            self.fail("extra files should of been deleted")

        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6_2"
        try:
            x3dna.generate_ref_frame(path)
            self.fail()
        except IOError:
            pass
        except:
            self.fail("got an unexpected error")

        # clean up
        os.remove("ref_frames.dat")

    def test_get_basepair_info(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6"
        x3dna = rnamake.x3dna.X3dna()
        basepairs = x3dna.get_basepair_info(path)
        if not os.path.isfile("ref_frames.dat"):
            self.fail("ref_frames.dat file should of been generated")
        if not os.path.isfile("p4p6_dssr.out"):
            self.fail("p4p6_dssr.out file should of been generated")

        os.remove("ref_frames.dat")
        os.remove("p4p6_dssr.out")

    def _test_get_basepair_info_compare(self):
        try:
            import redesign.motif
        except:
            self.skipTest("cannot load old REDESIGN package")

        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6"
        x3dna = rnamake.x3dna.X3dna()
        basepairs = x3dna.get_basepair_info(path)

    def test_generate_dssr_file(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6"
        x3dna = rnamake.x3dna.X3dna()
        x3dna.generate_dssr_file(path)
        if not os.path.isfile("p4p6_dssr.out"):
            self.fail("did not create _dssr.out file")
        os.remove("p4p6_dssr.out")

    def test_x3dna_residues(self):
        res1 = rnamake.x3dna.Residue(103,"A")
        res2 = rnamake.x3dna.Residue(103,"A")
        if res1 != res2:
            self.fail("residue equality is not working")

        res3 = rnamake.x3dna.Residue(102,"A")
        if res1 == res3:
            self.fail("residues should not be equal")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
