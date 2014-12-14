import unittest
import os
import rnamake.motif
import rnamake.settings
import rnamake.motif_library
import rnamake.motif_type
import rnamake.transform
import util
import numerical
import numpy as np

class MotifUnittest(unittest.TestCase):

    def setUp(self):
        path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
        self.motif = util.supress_log_output(rnamake.motif.Motif, path)

    def test_creation(self):

        m = rnamake.motif.Motif()

    def test_creation_mdir(self):
        try:
            path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
            m =  util.supress_log_output(rnamake.motif.Motif, path)
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
        if bp != found[0]:
            self.fail("did not retreive correct basepair")

        found = m.get_basepair(uuid1=bp.res1.uuid, uuid2=bp.res2.uuid)
        if bp != found[0]:
            self.fail("did not retreive correct basepair")

    def test_get_beads(self):
        m = self.motif
        beads = m.get_beads()
        org_count = len(beads)

        beads = m.get_beads(excluded_res=[m.structure.residues()[1]])
        diff = org_count - len(beads)
        if diff != 3:
            self.fail("did not exclude res properly")

        beads = m.get_beads([m.basepairs[0]])
        diff = org_count - len(beads)
        if diff != 5:
            self.fail("did not exclude ends properly")

    def test_to_str(self):
        m = self.motif
        s = m.to_str()
        m1 = rnamake.motif.str_to_motif(s)
        if len(m1.residues()) != 157:
            self.fail("did not copy all residues correctly")

    def test_copy(self):
        m = self.motif
        cm = m.copy()

    def test_align(self):
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library.MotifLibrary(mtype)
        m1 = mlib.get_motif("HELIX.IDEAL")
        m2 = mlib.get_motif("HELIX.IDEAL")
        #print m1.ends[1].state().r

        rnamake.motif.align_motif(m1.ends[1], m2.ends[0], m2)
        #m1.to_pdb("m1.pdb")
        #m2.to_pdb("m2.pdb")

    def test_transform(self):
        path = "/Users/josephyesselman/projects/REDESIGN/redesign/tests/p4p6"
        m = rnamake.motif.Motif(path)
        r = np.random.random([3,3])
        d = np.random.random([3])
        t = rnamake.transform.Transform(r, d)
        old_r = m.basepairs[0].state().r
        m.transform(t)
        new_r = m.basepairs[0].state().r
        if numerical.are_matrices_equal(old_r, new_r):
            self.fail("rotations should be different")

        m.reset()
        new_r = m.basepairs[0].state().r
        if not numerical.are_matrices_equal(old_r, new_r):
            self.fail("rotations should be different")










def main():
    unittest.main()

if __name__ == '__main__':
    main()
