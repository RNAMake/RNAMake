import unittest
import rnamake.motif_scorer
import rnamake.resource_manager
import redesign.motif_library

class MotifScorerUnittest(unittest.TestCase):

    def test_creation(self):
        scorer = rnamake.motif_scorer.MotifScorer()

    def test_score(self):
        rm =  rnamake.resource_manager.ResourceManager()
        scorer = rnamake.motif_scorer.MotifScorer()
        mlib = redesign.motif_library.MotifLibrary()

        score = scorer.score(rm.get_motif("HELIX.IDEAL"))
        score2 = mlib.get_motif("HELIX.IDEAL").score
        print score, score2

def main():
    unittest.main()

if __name__ == '__main__':
    main()
