import unittest
import os
import rnamake.setup.motif_library as motif_library
import rnamake.motif_type
import rnamake.settings

class MotifLibraryUnittest(unittest.TestCase):

    def setUp(self):
        if not os.path.isdir(rnamake.settings.MOTIF_DIRS + "/helices/"):
            self.skipTest("do not have the motif pdb files")

    def test_creation(self):
        mtype = rnamake.motif_type.TWOWAY
        mlib = motif_library.MotifLibrary(mtype)
        if len(mlib.motif_paths) == 0:
            self.fail("did not load motif paths properly")

        try:
            mlib = motif_library.MotifLibrary()
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("did not get the error I expected")

        path = rnamake.settings.MOTIF_DIRS + "two_ways/unique_7.dat"
        mlib = motif_library.MotifLibrary(libfile=path)
        if (mlib.motif_paths) == 0:
            self.fail("did not load motifs paths properly")

    def test_get_motif(self):
        mtype = rnamake.motif_type.HELIX
        mlib = motif_library.MotifLibrary(mtype)
        m = mlib.get_motif("HELIX.IDEAL")

        try:
            mlib.get_motif("HELIX")
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("did not get the error I expected")

    def test_get_motif_2(self):
        mtype = rnamake.motif_type.HELIX
        mlib = motif_library.MotifLibrary(mtype)
        m = mlib.get_motif("HELIX.IDEAL")
        m1 = mlib.get_motif("HELIX.IDEAL")

    def test_score(self):
        mtype = rnamake.motif_type.HELIX
        mlib = motif_library.MotifLibrary(mtype)
        m = mlib.get_motif("HELIX.IDEAL")
        if m.score != -3.8:
            self.fail("did not get the correct score")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
