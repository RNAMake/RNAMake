import unittest

from rnamake import sqlite_library

class MotifEnsembleUnittest(unittest.TestCase):

    def setUp(self):
        self.mlib = sqlite_library.MotifSqliteLibrary("twoways")

    def test_creation(self):
        motifs = []
        energies = []

        







def main():
    unittest.main()

if __name__ == '__main__':
    main()
