import unittest
import numpy as np
import rnamake.atom


class AtomUnittest(unittest.TestCase):

    def test_creation(self):
        atom = rnamake.atom.Atom("H1", np.array([0, 1, 2]))

    def test_slots(self):
        atom = rnamake.atom.Atom("H1", np.array([0, 1, 2]))
        self.assertRaises(AttributeError,atom.v1 = 1)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
