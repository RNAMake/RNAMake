import unittest
import rnamake.ss_tree as ss_tree

class SSTreeUnittest(unittest.TestCase):

    def test_creation(self):
        sstree = ss_tree.SS_Tree("((..(+)))", "GGUUC+GCC")
        #sstree = ss_tree.SS_Tree("(.(+).)", "GUG+CUC")

        #seq = "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA"
        #ss  = "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...))))).."
        #sstree = ss_tree.SS_Tree("((...(.).(.)...))", "GGAAAGACAGACAAACC")
        #sstree = ss_tree.SS_Tree(ss, seq)
        #for n in sstree:
        #    print n.data.ss_data

    def test_seq_from_nodes(self):
        sstree = ss_tree.SS_Tree("(.(+).)", "GUG+CUC")
        nodes = sstree.tree.nodes[0:2]
        ss_data = sstree.seq_from_nodes(nodes)
        print ss_data



def main():
    unittest.main()

if __name__ == '__main__':
    main()
