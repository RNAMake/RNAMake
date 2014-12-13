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

def main():
    unittest.main()

if __name__ == '__main__':
    main()
