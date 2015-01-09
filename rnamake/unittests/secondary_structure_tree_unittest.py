import unittest
import rnamake.secondary_structure_tree

class SecondaryStructureTreeUnittest(unittest.TestCase):

    def test_creation(self):
        ss  = ".....((..((.(......)))....(((((((....).)))..))).(((....).))..)).....(((((((....)))))))....................."
        seq = "GGAAAGCAAGGACGAAUAAGCCAUAACCAGAGCGAAAGACUCAAUGGAGCCGAAAGAGCAAGCAAUAACUGAUGCUUCGGCAUCAGAAAAGAAACAACAACAACAAC"
        ss_tree = rnamake.secondary_structure_tree.SecondaryStructureTree(ss,seq)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
