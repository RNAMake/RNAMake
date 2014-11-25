import unittest
import numpy as np
import rnamake.atom


class AtomUnittest(unittest.TestCase):

    def test_creation(self):
        try:
            rnamake.atom.Atom("H1", np.array([0, 1, 2]))
        except:
            self.fail("failed to initialize object")

    def test_slots(self):
        atom = rnamake.atom.Atom("H1", np.array([0, 1, 2]))

        try:
            atom.v1 = 1
        except AttributeError:
            pass
        except:
            self.fail("did not expect this error")

    def test_to_str(self):
        atom = rnamake.atom.Atom("H1", np.array([0, 1, 2]))
        string = atom.to_str()
        if string != "H1 0 1 2":
            print string
            self.fail("did not get correct string")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
