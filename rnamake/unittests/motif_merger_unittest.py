import unittest
from rnamake import sqlite_library, motif, motif_merger
import warnings

class MotifMergerUnittests(unittest.TestCase):

    def setUp(self):
        self.mlib_twoway   = sqlite_library.MotifSqliteLibrary("twoway")
        self.mlib_helix    = sqlite_library.MotifSqliteLibrary("ideal_helices")
        self.mlib_bp_steps = sqlite_library.MotifSqliteLibrary("new_bp_steps")

    def test_simple(self):
        m1 = self.mlib_helix.get(name="HELIX.IDEAL")
        m2 = self.mlib_helix.get(name="HELIX.IDEAL")

        motif.align_motif(m1.ends[1].state(), m2.ends[0], m2)
        mm = motif_merger.MotifMerger()
        mm.add_motif(m1)
        mm.add_motif(m2, m2.ends[0], m1, m1.ends[1])

        rna_struc = mm.get_structure()
        self.failIf(rna_struc.get_residue(uuid=m1.ends[0].res1.uuid) is None)
        self.failUnless(rna_struc.get_residue(uuid=m1.ends[1].res1.uuid) is None)

        ss = mm.secondary_structure()
        self.failUnless(ss.sequence() == "CCC&GGG")
        self.failUnless(ss.dot_bracket() == "(((&)))")

    def test_conserve_sequence_indentity_with_twoway(self):
        m1 = self.mlib_helix.get(name="HELIX.IDEAL")
        m2 = self.mlib_twoway.get(name="TWOWAY.1GID.12")

        motif.align_motif(m1.ends[1].state(), m2.ends[0], m2)
        mm = motif_merger.MotifMerger()
        mm.add_motif(m1)
        mm.add_motif(m2, m2.ends[0], m1, m1.ends[1])

        rna_struc = mm.get_structure()
        self.failIf(rna_struc.get_residue(uuid=m1.ends[0].res1.uuid) is None)
        self.failUnless(rna_struc.get_residue(uuid=m1.ends[1].res1.uuid) is None)

        m1 = self.mlib_helix.get(name="HELIX.IDEAL")
        m2 = self.mlib_twoway.get(name="TWOWAY.1GID.12")

        motif.align_motif(m2.ends[1].state(), m1.ends[0], m1)
        mm = motif_merger.MotifMerger()
        mm.add_motif(m2)
        mm.add_motif(m1, m1.ends[0], m2, m2.ends[1])
        rna_struc = mm.get_structure()

        # ideal helices sequence gets overwritten by twoway junction
        self.failIf(rna_struc.get_residue(uuid=m1.ends[0].res1.uuid) is not None)
        # should instead contain twoway junction basepair
        self.failIf(rna_struc.get_residue(uuid=m2.ends[1].res1.uuid) is None)

    def test_sequence_indenity_conflict(self):
        warnings.simplefilter("default")

        m1 = self.mlib_bp_steps.get(name="BP.1")
        m2 = self.mlib_bp_steps.get(name="BP.2")

        motif.align_motif(m1.ends[1].state(), m2.ends[0], m2)
        mm = motif_merger.MotifMerger()
        mm.add_motif(m1)

        # does not run in ./run_unittests
        with warnings.catch_warnings(record=True) as w:
            mm.add_motif(m2, m2.ends[0], m1, m1.ends[1])
            self.failUnless(len(w) > 0)

        m1 = self.mlib_bp_steps.get(end_id="GG_LL_CC_RR")
        m2 = self.mlib_bp_steps.get(end_id="GA_LL_UC_RR")
        motif.align_motif(m1.ends[1].state(), m2.ends[0], m2)
        mm = motif_merger.MotifMerger()
        mm.add_motif(m1)
        mm.add_motif(m2, m2.ends[0], m1, m1.ends[1])

        warnings.simplefilter("ignore")



def main():
    unittest.main()

if __name__ == '__main__':
    main()
