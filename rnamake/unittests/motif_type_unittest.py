import unittest
import rnamake.motif_type

class MotifTypeUnittest(unittest.TestCase):

    def test_type_to_str(self):
        mtype = rnamake.motif_type.TWOWAY
        s = rnamake.motif_type.type_to_str(mtype)
        if s != "TWOWAY":
            self.fail("did not get correct string")

    def test_str_to_type(self):
        s = "TWOWAY"
        mtype = rnamake.motif_type.str_to_type(s)
        if mtype != rnamake.motif_type.TWOWAY:
            self.fail("did not get correct type")

def main():
    unittest.main()

if __name__ == '__main__':
    main()
