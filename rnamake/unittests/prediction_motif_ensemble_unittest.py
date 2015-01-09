import unittest
import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree

class MotifEnsembleUnittest(unittest.TestCase):

    def test_creation(self):
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        me2 = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        met = motif_ensemble_tree.MotifEnsembleTree(me)
        met.add_ensemble(me2)



def main():
    unittest.main()

if __name__ == '__main__':
    main()
