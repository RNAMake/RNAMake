import unittest
import rnamake.motif_library
import rnamake.motif_type
import rnamake.settings

class MotifLibraryUnittest(unittest.TestCase):

    def test_creation(self):
        mtype = rnamake.motif_type.TWOWAY
        mlib = rnamake.motif_library.MotifLibrary(mtype)
        if len(mlib.motif_paths) == 0:
            self.fail("did not load motif paths properly")

        try:
            mlib = rnamake.motif_library.MotifLibrary()
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("did not get the error I expected")

        path = rnamake.settings.MOTIF_DIRS + "two_ways/unique_7.dat"
        mlib = rnamake.motif_library.MotifLibrary(libfile=path)
        if (mlib.motif_paths) == 0:
            self.fail("did not load motifs paths properly")


    def test_get_motif(self):
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library.MotifLibrary(mtype)
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
        mlib = rnamake.motif_library.MotifLibrary(mtype)
        m = mlib.get_motif("HELIX.IDEAL")
        m1 = mlib.get_motif("HELIX.IDEAL")
        print m.ends[0].d()
        print m1.ends[0].d()



def main():
    unittest.main()

if __name__ == '__main__':
    main()
