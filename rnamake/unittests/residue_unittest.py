import unittest
import rnamake.residue_type
import rnamake.residue


class ResidueUnittest(unittest.TestCase):

    def test_creation(self):
        gtype = rnamake.residue_type.get_rtype("GUA")
        try:
            res = rnamake.residue.Residue(gtype, "GUA", 1, "A")
        except:
            self.fail("cannot creat Residue object sucessfully")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
