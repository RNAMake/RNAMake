import unittest
import rnamake.prediction.ensemble_predictor as ensemble_predictor
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree
import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.secondary_structure_tree as secondary_structure_tree
import rnamake.settings as settings
import rnamake.motif as motif
import random

class EnsemblePredictorUnittest(unittest.TestCase):

    def test_creation(self):
        return
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        me2 = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        met = motif_ensemble_tree.MotifEnsembleTree(me)
        met.add_ensemble(me2)
        ep = ensemble_predictor.EnsemblePredictor(met=met)

    def test_creation_ss_structure(self):
        ss  = "(((())))"
        seq = "ACGCGCGU"
        ep = ensemble_predictor.EnsemblePredictor(seq, ss)
        ep.sample(output_pdbs=0)

    def test_creation_ss_structure_2(self):
        seq, ss = self._generate_helix(size=100)
        #print seq
        #print ss
        ep = ensemble_predictor.EnsemblePredictor(seq, ss)
        ep.sample(output_pdbs=0, steps=50)

    def _generate_helix(self, size=100):
        seq = ""
        ss = ""

        steps = ["AU", "UA", "CG", "CG"]
        for i in range(size):
            step = random.choice(steps)
            seq = step[0] + seq  + step[1]
            ss = "(" + ss + ")"
        return seq, ss


    def test_tecto_rnas(self):
        return
        seq1 = "GGACUAGGAUAUGGAAGAUCCUCGGGAACGAGGAUCUUCCUAAGUCCUAG"
        ss1  = "...(((((((..((((((((((((....))))))))))))...)))))))"
        seq2 = "CUAGGAAUCUGGAAGAUCCUCGGAAACGAGGAUCUUCCUGUGUCCUAG"
        ss2  = "((((((....((((((((((((....))))))))))))....))))))"
        seq = seq1+seq2
        ss = ss1+ss2
        ss_tree = secondary_structure_tree.SecondaryStructureTree(ss, seq)
        #print ss_tree.nodes[0].children
        motifs = []
        motifs.append(motif.Motif(settings.UNITTEST_PATH + "/resources/GAAA_tetraloop"))
        motifs.append(motif.Motif(settings.UNITTEST_PATH + "/resources/GGAA_tetraloop"))
        ep = ensemble_predictor.EnsemblePredictor(seq, ss, motifs)



def main():
    unittest.main()

if __name__ == '__main__':
    main()
