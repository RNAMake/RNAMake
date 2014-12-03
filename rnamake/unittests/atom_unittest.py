import unittest
import numpy as np
import rnamake.atom
import rnamake.io


class AtomUnittest(unittest.TestCase):

    def test_creation(self):
        """
        tests basic initiation of the atom class
        """
        try:
            rnamake.atom.Atom("H1", np.array([0, 1, 2]))
        except:
            self.fail("failed to initialize object")

    def test_slots(self):
        """
        tests to make sure that no other attributes can be added to atoms
        if they are an AttributeError is thrown
        """
        atom = rnamake.atom.Atom("H1", np.array([0, 1, 2]))

        try:
            atom.v1 = 1
        except AttributeError:
            pass
        except:
            self.fail("did not expect this error")

    def test_to_str(self):
        """
        tests whether to_str() formats data correctly
        """
        atom = rnamake.atom.Atom("H1", np.array([0, 1, 2]))
        string = atom.to_str()
        if string != "H1 0 1 2":
            print string
            self.fail("did not get correct string")

    def test_to_pdb_str(self):
        atom = rnamake.atom.Atom("H1", np.array([1, 2, 3]))
        string = atom.to_pdb_str()
        refs = "ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n"
        if string != refs:
            print
            print "Actual String  ", string,
            print "Expected String", refs
            self.fail("did not get the correct pdb string")
        string = atom.to_pdb_str(10)
        refs = "ATOM     10  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n"
        if string != refs:
            print
            print "Actual String  ", string,
            print "Expected String", refs
            self.fail("did not get the correct pdb string")

    def test_str_to_atom(self):
        atom = rnamake.atom.Atom("H1", np.array([1.0, 2.0, 3.0]))
        string = atom.to_str()
        atom2 = rnamake.io.str_to_atom(string)
        self.assertListEqual(atom.coords.tolist(), atom2.coords.tolist())
        self.assertEqual(atom.name, atom2.name)

def main():
    unittest.main()

if __name__ == '__main__':
    main()
