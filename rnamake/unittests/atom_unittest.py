import unittest
import numpy as np

from rnamake import util, exceptions, all_atom
from instances import transform_instances
import is_equal


class AtomUnittest(unittest.TestCase):

    def setUp(self):
        self.a = all_atom.Atom("H1", np.array([0, 1, 2]))

    def test_creation(self):
        """
        tests basic initiation of the atom class
        """
        try:
            all_atom.Atom("H1", np.array([0, 1, 2]))
        except:
            self.fail("failed to initialize object")

        with self.assertRaises(exceptions.AtomException):
            all_atom.Atom("H1", [0, 1, 2])

        a = self.a

        # cannot override internal value if using the @getter
        coords = a.get_coords()
        coords[0] = 10
        self.failUnless(a.get_coords()[0] == 0)

    def test_slots(self):
        """
        tests to make sure that no other attributes can be added to atoms
        if they are an AttributeError is thrown
        """
        a = self.a

        try:
            a.v1 = 1
        except AttributeError:
            pass
        except:
            self.fail("did not expect this error")

    def test_copy(self):
        a = self.a
        copy_a = all_atom.Atom.copy(a)

        copy_a.move(np.array([1,1,1]))
        diff = util.distance(a.get_coords(),copy_a.get_coords())
        if diff < 0.1:
            self.fail()

        copy_a = a.get_copy()
        copy_a.move(np.array([1, 1, 1]))
        diff = util.distance(a.get_coords(), copy_a.get_coords())
        if diff < 0.1:
            self.fail()

    def test_to_str(self):
        """
        tests whether to_str() formats data correctly
        """
        a = self.a
        string = a.get_str()
        if string != "H1 0.0 1.0 2.0":
            print string
            self.fail("did not get correct string")

    def test_to_pdb_str(self):
        a = all_atom.Atom("H1", np.array([1, 2, 3]))
        string = a.get_pdb_str()
        refs = "ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n"
        if string != refs:
            print
            print "Actual String  ", string,
            print "Expected String", refs
            self.fail("did not get the correct pdb string")
        string = a.get_pdb_str(10)
        refs = "ATOM     10  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n"
        if string != refs:
            print
            print "Actual String  ", string,
            print "Expected String", refs
            self.fail("did not get the correct pdb string")

    def test_str_to_atom(self):
        a = self.a
        s = a.get_str()
        a2 = all_atom.Atom.from_str(s)
        self.assertListEqual(a.get_coords().tolist(), a2.get_coords().tolist())
        self.assertEqual(a.get_name(), a2.get_name())

        s = "N"
        with self.assertRaises(exceptions.AtomException):
            all_atom.Atom.from_str(s)

    def test_move(self):
        a = self.a
        a.move(np.array([0, 0, 1]))

        new_a = all_atom.Atom("H1", np.array([0, 1, 3]))
        self.failUnless(is_equal.are_atom_equal(a, new_a))

    def test_transform(self):
        a = all_atom.Atom.copy(self.a)
        t = transform_instances.transform_indentity()
        a.transform(t)

        self.failUnless(is_equal.are_atom_equal(a, self.a))

        t = transform_instances.transform_random()
        a.transform(t)

        self.failUnless(not is_equal.are_atom_equal(a, self.a))





def main():
    unittest.main()

if __name__ == '__main__':
    main()
