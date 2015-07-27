import unittest
import rnamake.motif_scorer
import rnamake.resource_manager

class MotifScorerUnittest(unittest.TestCase):

    def test_creation(self):
        scorer = rnamake.motif_scorer.MotifScorer()

    def test_score(self):
        rm =  rnamake.resource_manager.ResourceManager()
        scorer = rnamake.motif_scorer.MotifScorer()
        score = scorer.score(rm.get_motif(name="HELIX.IDEAL"))
        if score != -3.8:
            self.fail("did not get the right score")

def main():
    unittest.main()

if __name__ == '__main__':
    main()
