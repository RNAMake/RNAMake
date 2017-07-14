import unittest
from rnamake import sqlite_library, motif, motif_merger, motif_factory
import warnings

class MotifMergerUnittests(unittest.TestCase):

    def setUp(self):
        self.mf = motif_factory.MotifFactory()
        self.mlib_twoway   = sqlite_library.MotifSqliteLibrary("twoway")
        self.mlib_helix    = sqlite_library.MotifSqliteLibrary("ideal_helices")
        self.mlib_bp_steps = sqlite_library.MotifSqliteLibrary("bp_steps")

    def test_simple(self):
        m1 = self.mlib_helix.get(name="HELIX.IDEAL")
        m2 = self.mlib_helix.get(name="HELIX.IDEAL")

        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        mm = motif_merger.MotifMerger(self.mf)
        mm.add_motif(m1)
        mm.add_motif(m2, m2.get_end(0), m1, m1.get_end(1))

        rna_struc = mm.get_merged_structure()

        self.failIf(rna_struc.get_residue(uuid=m1.get_end(0).res1_uuid) is None)
        self.failUnless(rna_struc.get_residue(uuid=m1.get_end(1).res1_uuid) is None)

        ss = mm.get_merged_secondary_structure()
        self.failUnless(ss.sequence() == "GGG&CCC")
        self.failUnless(ss.dot_bracket() == "(((&)))")

    def test_conserve_sequence_indentity_with_twoway(self):
        m1 = self.mlib_helix.get(name="HELIX.IDEAL")
        m2 = self.mlib_twoway.get(name="TWOWAY.1GID.12")

        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        mm = motif_merger.MotifMerger(self.mf)
        mm.add_motif(m1)
        mm.add_motif(m2, m2.get_end(0), m1, m1.get_end(1))

        rna_struc = mm.get_merged_structure()
        self.failIf(rna_struc.get_residue(uuid=m1.get_end(0).res1_uuid) is None)
        self.failUnless(rna_struc.get_residue(uuid=m1.get_end(1).res1_uuid) is None)

        m1 = self.mlib_helix.get(name="HELIX.IDEAL")
        m2 = self.mlib_twoway.get(name="TWOWAY.1GID.12")

        motif.align_motif(m2.get_end(1), m1.get_end(0), m1)
        mm = motif_merger.MotifMerger(self.mf)
        mm.add_motif(m2)
        mm.add_motif(m1, m1.get_end(0), m2, m2.get_end(1))
        rna_struc = mm.get_merged_structure()

        # ideal helices sequence gets overwritten by twoway junction
        self.failIf(rna_struc.get_residue(uuid=m1.get_end(0).res1_uuid) is not None)
        # should instead contain twoway junction basepair
        self.failIf(rna_struc.get_residue(uuid=m2.get_end(1).res1_uuid) is None)

    def test_sequence_indenity_conflict(self):
        warnings.simplefilter("default")

        m1 = self.mlib_bp_steps.get(end_id="GG_LL_CC_RR")
        m2 = self.mlib_bp_steps.get(end_id="AA_LL_UU_RR")

        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        mm = motif_merger.MotifMerger(self.mf)
        mm.add_motif(m1)

        # does not run in ./run_unittests
        with warnings.catch_warnings(record=True) as w:
            mm.add_motif(m2, m2.get_end(0), m1, m1.get_end(1))
            self.failUnless(len(w) > 0)

        m1 = self.mlib_bp_steps.get(end_id="GG_LL_CC_RR")
        m2 = self.mlib_bp_steps.get(end_id="GA_LL_UC_RR")
        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        mm = motif_merger.MotifMerger(self.mf)
        mm.add_motif(m1)
        mm.add_motif(m2, m2.get_end(0), m1, m1.get_end(1))

        warnings.simplefilter("ignore")

class MotifStateMergerUnittests(unittest.TestCase):

    def setUp(self):
        self.mf = motif_factory.MotifFactory()
        self.mlib_twoway   = sqlite_library.MotifSqliteLibrary("twoway")
        self.mlib_helix    = sqlite_library.MotifSqliteLibrary("ideal_helices")
        self.mlib_bp_steps = sqlite_library.MotifSqliteLibrary("bp_steps")

    def test_simple(self):
        m1 = self.mlib_helix.get(name="HELIX.IDEAL").get_state()
        m2 = self.mlib_helix.get(name="HELIX.IDEAL").get_state()

        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        msm = motif_merger.MotifStateMerger()
        msm.add_motif(m1)
        msm.add_motif(m2, m2.get_end(0), m1, m1.get_end(1))

        rna_struc = msm.get_merged_structure()

        self.failIf(rna_struc.get_residue(uuid=m1.get_end(0).res1_uuid) is None)
        self.failUnless(rna_struc.get_residue(uuid=m1.get_end(1).res1_uuid) is None)

        ss = msm.get_merged_secondary_structure()
        self.failUnless(ss.sequence() == "GGG&CCC")
        self.failUnless(ss.dot_bracket() == "(((&)))")

    def test_conserve_sequence_indentity_with_twoway(self):
        m1 = self.mlib_helix.get(name="HELIX.IDEAL").get_state()
        m2 = self.mlib_twoway.get(name="TWOWAY.1GID.12").get_state()

        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        msm = motif_merger.MotifStateMerger()
        msm.add_motif(m1)
        msm.add_motif(m2, m2.get_end(0), m1, m1.get_end(1))

        rna_struc = msm.get_merged_structure()
        self.failIf(rna_struc.get_residue(uuid=m1.get_end(0).res1_uuid) is None)
        self.failUnless(rna_struc.get_residue(uuid=m1.get_end(1).res1_uuid) is None)

        m1 = self.mlib_helix.get(name="HELIX.IDEAL").get_state()
        m2 = self.mlib_twoway.get(name="TWOWAY.1GID.12").get_state()

        motif.align_motif(m2.get_end(1), m1.get_end(0), m1)
        mm = motif_merger.MotifStateMerger()
        msm.add_motif(m2)
        msm.add_motif(m1, m1.get_end(0), m2, m2.get_end(1))
        rna_struc = msm.get_merged_structure()

        # ideal helices sequence gets overwritten by twoway junction
        self.failIf(rna_struc.get_residue(uuid=m1.get_end(0).res1_uuid) is not None)
        # should instead contain twoway junction basepair
        self.failIf(rna_struc.get_residue(uuid=m2.get_end(1).res1_uuid) is None)

    def _test_sequence_indenity_conflict(self):
        warnings.simplefilter("default")

        m1 = self.mlib_bp_steps.get(end_id="GG_LL_CC_RR")
        m2 = self.mlib_bp_steps.get(end_id="AA_LL_UU_RR")

        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        mm = motif_merger.MotifMerger(self.mf)
        mm.add_motif(m1)

        # does not run in ./run_unittests
        with warnings.catch_warnings(record=True) as w:
            mm.add_motif(m2, m2.get_end(0), m1, m1.get_end(1))
            self.failUnless(len(w) > 0)

        m1 = self.mlib_bp_steps.get(end_id="GG_LL_CC_RR")
        m2 = self.mlib_bp_steps.get(end_id="GA_LL_UC_RR")
        motif.align_motif(m1.get_end(1), m2.get_end(0), m2)
        mm = motif_merger.MotifMerger(self.mf)
        mm.add_motif(m1)
        mm.add_motif(m2, m2.get_end(0), m1, m1.get_end(1))

        warnings.simplefilter("ignore")



def main():
    unittest.main()

if __name__ == '__main__':
    main()
