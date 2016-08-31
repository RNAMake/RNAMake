import unittest

from rnamake import sqlite_library, settings, motif, motif_merger, motif_factory


class MotifMergerTwowayTests(unittest.TestCase):

    def setUp(self):
        self.mlib = sqlite_library.MotifSqliteLibrary("twoway")
        self.mlib.load_all()

        path = settings.RESOURCES_PATH + "/motifs/base.motif"
        self.base_motif = motif.file_to_motif(path)

    def test_merge(self):

        mm = motif_merger.MotifMerger()
        mm.add_motif(self.base_motif)
        m2 = self.base_motif.copy()
        m2.new_res_uuids()
        for m in self.mlib.all():
            mm.add_motif(m, m.ends[0], self.base_motif, self.base_motif.ends[1])
            s = mm.get_structure()

            self.failUnless(len(s.chains()) == 2)

            mm.remove_motif(m)

class MotifMergerNwayTests(unittest.TestCase):

    def setUp(self):
        self.mlib = sqlite_library.MotifSqliteLibrary("nway")
        self.mlib.load_all()

        path = settings.RESOURCES_PATH + "/motifs/base.motif"
        self.base_motif = motif.file_to_motif(path)

    def test_merge(self):

        mm = motif_merger.MotifMerger()
        mm.add_motif(self.base_motif)
        m2 = self.base_motif.copy()
        m2.new_res_uuids()
        for m in self.mlib.all():

            mm.add_motif(m, m.ends[0], self.base_motif, self.base_motif.ends[1])
            s = mm.get_structure()

            self.failUnless(len(s.chains()) == len(m.chains()))

            mm.remove_motif(m)

def main():
    unittest.main()

if __name__ == '__main__':
    main()