import unittest
import rnamake.motif_tree_merger
import instance

class MotifTreeMergerUnittest(unittest.TestCase):

    def test_creation(self):
        mt = instance.simple_mt()
        merger = rnamake.motif_tree_merger.MotifTreeMerger()
        mt.write_pdbs()
        motif = merger.merge(mt)
        motif.to_pdb()

def main():
    unittest.main()

if __name__ == '__main__':
    main()
