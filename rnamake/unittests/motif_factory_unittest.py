import unittest
import os
from rnamake import motif_factory, settings


class MotifFactoryUnittest(unittest.TestCase):
    def setUp(self):
        self.mf = motif_factory.MotifFactory()

    def test_creation(self):
        pass

    def test_motifs_from_file(self):
        path = settings.MOTIF_DIRS + "junctions/NWAY.1GID.0"
        motifs = self.mf.motifs_from_file(path)
        self.failUnless(len(motifs) == 3)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
