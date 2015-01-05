import unittest
import rnamake.prediction.setup
import rnamake.motif_type
import rnamake.motif_library
import util

class SetupUnittest(unittest.TestCase):

    def test_get_bp_setup(self):
        target = ["GC","GC"]
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library.MotifLibrary(mtype)
        mlib.load_all()

        rnamake.prediction.setup.get_bp_set_prediction_lib(target, mlib)

def main():
    unittest.main()

if __name__ == '__main__':
    main()
