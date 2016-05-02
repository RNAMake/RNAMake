import unittest

from rnamake import exceptions, motif_type

class MotifTypeUnittest(unittest.TestCase):

    def test_type_to_str(self):
        mtype = motif_type.TWOWAY
        s = motif_type.type_to_str(mtype)
        if s != "TWOWAY":
            self.fail("did not get correct string")

        with self.assertRaises(exceptions.MotifTypeException):
            motif_type.type_to_str("FAKE")

    def test_str_to_type(self):
        s = "TWOWAY"
        mtype = motif_type.str_to_type(s)
        if mtype != motif_type.TWOWAY:
            self.fail("did not get correct type")

def main():
    unittest.main()

if __name__ == '__main__':
    main()
