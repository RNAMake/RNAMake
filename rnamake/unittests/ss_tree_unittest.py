import unittest
import rnamake.ss_tree as ss_tree

class SSTreeUnittest(unittest.TestCase):

    def test_creation(self):
        #sstree = ss_tree.SS_Tree("((+))", "GG+CC")


        #sstree = ss_tree.SS_Tree("(.(+).)", "GUG+CUC")

        seq = "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA"
        ss  = "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...))))).."
        #sstree = ss_tree.SS_Tree("((...(.).(.)...))", "GGAAAGACAGACAAACC")
        sstree = ss_tree.SS_Tree(ss, seq)

        for n in sstree:
            print n.data.seqs



def main():
    unittest.main()

if __name__ == '__main__':
    main()
