import unittest
import os
import shutil
import rnamake.x3dna

class X3dnaUnittest(unittest.TestCase):

    def test_generate_ref_frame(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb"
        x3dna = rnamake.x3dna.X3dna()
        x3dna.generate_ref_frame(path)
        if not os.path.isfile("ref_frames.dat"):
            self.fail("ref_frames.dat file should of been generated")
        if os.path.isfile("bestpairs.pdb"):
            self.fail("extra files should of been deleted")

        path = rnamake.settings.UNITTEST_PATH + "resources/motifs/p4p6_2"
        try:
            x3dna.generate_ref_frame(path)
            self.fail()
        except IOError:
            pass
        except:
            self.fail("got an unexpected error")

        # clean up
        if not os.path.isfile("ref_frames.dat"):
            self.fail("did not produce ref_frames.dat")

        os.remove("ref_frames.dat")

    def test_generate_dssr_file(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb"
        x3dna = rnamake.x3dna.X3dna()
        x3dna.generate_dssr_file(path)
        if not os.path.isfile("p4p6_dssr.out"):
            self.fail("did not create _dssr.out file")

        # clean up
        if not os.path.isfile("p4p6_dssr.out"):
            self.fail("did not produce p4p6_dssr.out")

        os.remove("p4p6_dssr.out")

    def test_get_basepairs(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb"
        x3dna = rnamake.x3dna.X3dna()
        basepairs = x3dna.get_basepairs(path)
        if len(basepairs) != 79:
            self.fail("did not produced the expected number of basepairs")


        #os.remove("ref_frames.dat")
        #os.remove("p4p6_dssr.out")

    def test_get_basepairs_compare(self):
        try:
            import redesign.motif
        except:
            self.skipTest("cannot load old REDESIGN package")

        path = rnamake.settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb"
        x3dna = rnamake.x3dna.X3dna()
        basepairs = x3dna.get_basepairs(path)

        try:
            path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
            m = redesign.motif.Motif(path)
        except:
            self.skipTest("something is wrong cannot load p4p6 with REDESIGN")

        for bp in m.base_pairs:
            found = 0
            for xbp in basepairs:
                if bp.res1.num == xbp.res1.num and bp.res2.num == xbp.res2.num:
                    found = 1
                    break
                if bp.res1.num == xbp.res2.num and bp.res2.num == xbp.res1.num:
                    found = 1
                    break

            if not found:
                print bp.res1,bp.res2
                self.fail()

    def test_x3dna_residues(self):
        res1 = rnamake.x3dna.Residue(103,"A")
        res2 = rnamake.x3dna.Residue(103,"A")
        if res1 != res2:
            self.fail("residue equality is not working")

        res3 = rnamake.x3dna.Residue(102,"A")
        if res1 == res3:
            self.fail("residues should not be equal")

    def test_get_motifs(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb"
        x3dna = rnamake.x3dna.X3dna()
        motifs = x3dna.get_motifs(path)

    def test_get_basepairs_local(self):
        path1 = rnamake.settings.UNITTEST_PATH + "resources/motifs/p4p6/p4p6.pdb"
        path = rnamake.settings.UNITTEST_PATH + "p4p6.pdb"
        shutil.copy(path1, path)
        x3dna = rnamake.x3dna.X3dna()
        basepairs = x3dna.get_basepairs(path)
        os.remove("ref_frames.dat")
        os.remove('p4p6_dssr.out')
        os.remove("p4p6.pdb")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
