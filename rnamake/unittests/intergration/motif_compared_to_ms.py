import unittest
from rnamake import sqlite_library, settings, motif, residue_type, util


class BasicLibrariesUnittests(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()

        path = settings.RESOURCES_PATH + "/motifs/base.motif"
        self.base_motif = motif.file_to_motif(path)
        self.added_motif = motif.file_to_motif(path)

        self.base_ms = self.base_motif.get_state()
        self.added_ms = self.added_motif.get_state()

    def _test_correct_build(self, mlib, ms_lib):

        for m in mlib.all():
            m1 = motif.get_aligned_motif(self.base_motif.ends[1], m.ends[0], m)
            m2 = motif.get_aligned_motif(m1.ends[1],
                                         self.added_motif.ends[0],
                                         self.added_motif)

            ms = ms_lib.get(name=m.name,
                            end_name=m.ends[0].name(),
                            end_id=m.end_ids[0])

            ms1 = motif.get_aligned_motif_state_single(self.base_ms.end_states[1], ms)
            ms2 = motif.get_aligned_motif_state_single(ms1.end_states[1], self.added_ms)

            #print m2.ends[1].d()
            #print self.added_ms.end_states[1].d
            diff = ms2.end_states[1].diff(m2.ends[1].state())
            if diff > 0.01:
                self.fail(m.name + " did not give the same answer as its state")



    def test_correct_build_twoway(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        mlib.load_all()
        ms_lib = sqlite_library.MotifStateSqliteLibrary("twoway")
        ms_lib.load_all()
        self._test_correct_build(mlib, ms_lib)

    def test_correct_build_nway(self):
        mlib = sqlite_library.MotifSqliteLibrary("nway")
        mlib.load_all()
        ms_lib = sqlite_library.MotifStateSqliteLibrary("nway")
        ms_lib.load_all()
        self._test_correct_build(mlib, ms_lib)




def main():
    unittest.main()

if __name__ == '__main__':
    main()