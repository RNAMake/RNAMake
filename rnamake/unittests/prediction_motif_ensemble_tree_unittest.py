import unittest
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree
import rnamake.prediction.motif_ensemble as motif_ensemble

class MotifEnsembleTreeUnittest(unittest.TestCase):

    def test_creation(self):
        met = motif_ensemble_tree.MotifEnsembleTree()

    def test_add(self):
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        met = motif_ensemble_tree.MotifEnsembleTree()
        met.add_ensemble(me)

    def test_get_mtst(self):
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        me2 = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        met = motif_ensemble_tree.MotifEnsembleTree()
        met.add_ensemble(me)
        met.add_ensemble(me2)
        mtst = met.get_mtst()
        print len(mtst.nodes)
        mtst.nodes_to_pdbs()

def main():
    unittest.main()

if __name__ == '__main__':
    main()
