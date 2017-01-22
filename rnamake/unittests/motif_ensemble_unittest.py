import unittest

from rnamake import sqlite_library, motif, motif_ensemble, residue_type
import is_equal

class MotifEnsembleUnittest(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()
        self.mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        m1 = self.mlib.get(name="HELIX.IDEAL")
        m2 = motif.Motif.copy(m1)

        motifs = [m1, m2]
        energies = [1, 1]

        self.me = motif_ensemble.MotifEnsemble(motifs, energies)

    def test_creation(self):
        me = self.me
        self.failUnless(len(me) == 2)

    def test_to_str(self):
        me = self.me
        s = me.to_str()
        me2 = motif_ensemble.MotifEnsemble.from_str(s, self.rts)

        self.failUnless(len(me) == len(me2))
        self.failUnless(is_equal.are_rna_strucs_equal(me.get_member(0).motif,
                                                      me2.get_member(0).motif,
                                                      check_uuid=0))

        self.failUnless(is_equal.are_rna_strucs_equal(me.get_member(1).motif,
                                                      me2.get_member(1).motif,
                                                      check_uuid=0))

    def test_copy(self):
        me = self.me
        me2 = motif_ensemble.MotifEnsemble.copy(me)

        self.failUnless(len(me) == len(me2))
        self.failUnless(is_equal.are_rna_strucs_equal(me.get_member(0).motif,
                                                      me2.get_member(0).motif,
                                                      check_uuid=0))

        self.failUnless(is_equal.are_rna_strucs_equal(me.get_member(1).motif,
                                                      me2.get_member(1).motif,
                                                      check_uuid=0))

    def test_get_state(self):
        me = self.me
        mse = me.get_state()

        self.failUnless(len(me) == len(mse))

def main():
    unittest.main()

if __name__ == '__main__':
    main()
