import unittest
import rnamake.motif_tree_precomputer

class MotifTreePrecomputerUnittest(unittest.TestCase):

    def test_creation(self):
        mtp = rnamake.motif_tree_precomputer.MotifTreePrecomputer()


def main():
    unittest.main()

if __name__ == '__main__':
    main()
