import unittest
import rnamake.motif_library_sqlite
import rnamake.motif_type

class MotifLibrarySqliteUnittest(unittest.TestCase):

    def test_creation(self):
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library_sqlite.MotifLibrarySqlite(mtype)

    def test_get_motif(self):
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library_sqlite.MotifLibrarySqlite(mtype)
        m = mlib.get_motif("HELIX.IDEAL")

        try:
            mlib.get_motif("HELIX")
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("did not get the error I expected")

    def test_load_all(self):
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library_sqlite.MotifLibrarySqlite(mtype)
        mlib.load_all(limit=10)
        for mname, m in mlib.mdict.iteritems():
            print mname, m




def main():
    unittest.main()

if __name__ == '__main__':
    main()
