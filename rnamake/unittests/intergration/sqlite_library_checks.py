import unittest

from rnamake import sqlite_library, settings, motif, residue_type


class IdealHelicesUnittests(unittest.TestCase):

    def setUp(self):
        self.mlib =  sqlite_library.MotifSqliteLibrary("ideal_helices")
        self.mlib.load_all()

        self.base_motif = self.mlib.get(name="HELIX.IDEAL.1")
        self.added_motif = motif.Motif.copy(self.base_motif)

    def test_correct_build(self):

        for m in self.mlib.all():

            m1 = motif.get_aligned_motif(self.base_motif.get_end(1), m.get_end(0), m)
            self.failIf(motif.clash_between_motifs(self.base_motif, m1))

            m2 = motif.get_aligned_motif(m1.get_end(1), self.added_motif.get_end(0),
                                         self.added_motif)

            self.failIf(motif.clash_between_motifs(m1, m2))
            self.failIf(motif.clash_between_motifs(self.base_motif, m2))


class IdealReverseHelicesUnittests(unittest.TestCase):

    def setUp(self):
        self.mlib =  sqlite_library.MotifSqliteLibrary("ideal_helices_reversed")
        self.mlib.load_all()
        self.rts = residue_type.ResidueTypeSet()

        path = settings.RESOURCES_PATH + "/motifs/base.motif"
        self.base_motif = motif.file_to_motif(path)
        self.added_motif = motif.file_to_motif(path)

    def test_correct_build(self):

        mg = motif_graph.MotifGraph()
        mg.add_motif(self.base_motif)
        mg.increase_level()

        for m in self.mlib.all():
            mg.add_motif(m)
            mg.add_motif(self.added_motif)
            if len(mg) != 3:
                self.fail(m.name + " did not build correctly")


            mg.remove_node_level()

    def test_secondary_structure(self):

        m = self.mlib.get(name="HELIX.IDEAL.3")

        if m.sequence() == self.base_motif.sequence():
            self.fail("sequence should be flipped but now")


class BasicLibrariesUnittests(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()
        mlib =  sqlite_library.MotifSqliteLibrary("ideal_helices")

        self.base_motif = mlib.get(name="HELIX.IDEAL.1")
        self.added_motif = motif.Motif.copy(self.base_motif)

    def _test_correct_build(self, mlib):

        for m in mlib.all():
            m1 = motif.get_aligned_motif(self.base_motif.get_end(1), m.get_end(0), m)

            self.failIf(motif.clash_between_motifs(self.base_motif, m1))

            for i in range(1, m1.num_ends()):
                m2 = motif.get_aligned_motif(m1.get_end(i), self.added_motif.get_end(0),
                                             self.added_motif)

                self.failIf(motif.clash_between_motifs(m1, m2))
                self.failIf(motif.clash_between_motifs(self.base_motif, m2))

    def test_correct_build_twoway(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        mlib.load_all()
        self._test_correct_build(mlib)

    def test_correct_build_nway(self):
        mlib = sqlite_library.MotifSqliteLibrary("nway")
        mlib.load_all()
        self._test_correct_build(mlib)

    def test_correct_build_tcontact(self):
        mlib = sqlite_library.MotifSqliteLibrary("tcontact")
        mlib.load_all()
        self._test_correct_build(mlib)

    def test_correct_build_hairpin(self):
        mlib = sqlite_library.MotifSqliteLibrary("hairpin")
        mlib.load_all()
        self._test_correct_build(mlib)

    def test_correct_build_bp_steps(self):
        mlib = sqlite_library.MotifSqliteLibrary("bp_steps")
        mlib.load_all()
        self._test_correct_build(mlib)


class BPStepsUnittests(unittest.TestCase):

    def setUp(self):

        path = settings.RESOURCES_PATH + "/motifs/base.motif"
        self.rts = residue_type.ResidueTypeSet()
        self.base_motif = motif.file_to_motif(path)

    def _test_correct_build(self, me):

        for mem in me.members:
            m = mem.motif
            m1 = motif.get_aligned_motif(self.base_motif.ends[1], m.ends[0], m)

            if motif.clash_between_motifs(self.base_motif, m1):
                self.fail("clash")

            for i in range(1,len(m1.ends)):
                m2 = motif.get_aligned_motif(m1.ends[i], self.base_motif.ends[0],
                                             self.base_motif)

                if motif.clash_between_motifs(m1, m2):
                    self.fail("clash")
                if motif.clash_between_motifs(self.base_motif, m2):
                    self.fail("clash")

    def test_correct_build_bps(self):
        me_lib = sqlite_library.MotifEnsembleSqliteLibrary("bp_steps")
        me_lib.load_all()

        for me in me_lib.all():
            self._test_correct_build(me)




def main():
    unittest.main()

if __name__ == '__main__':
    main()