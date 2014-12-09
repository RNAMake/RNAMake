import unittest
import os
import rnamake.x3dna

class X3dnaUnittest(unittest.TestCase):

    def test_generate_ref_frame(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6"
        x3dna = rnamake.x3dna.X3dna()
        x3dna.generate_ref_frame(path)
        if not os.path.isfile("ref_frames.dat"):
            self.fail("ref_frames.dat file should of been generated")
        if os.path.isfile("bestpairs.pdb"):
            self.fail("extra files should of been deleted")

        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6_2"
        try:
            x3dna.generate_ref_frame(path)
            self.fail()
        except IOError:
            pass
        except:
            self.fail("got an unexpected error")






def main():
    unittest.main()

if __name__ == '__main__':
    main()
