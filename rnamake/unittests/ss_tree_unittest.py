import unittest
import rnamake.ss_tree as ss_tree

class SSTreeUnittest(unittest.TestCase):

    def test_creation(self):
        sstree = ss_tree.SS_Tree("GAG&CAG&CAC", "(.(&).(&).)")
        #sstree = ss_tree.SS_Tree("GGUUC+GCC", "((..(+)))")


        #sstree = ss_tree.SS_Tree("(.(+).)", "GUG+CUC")

        seq = "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA"
        ss  = "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...))))).."
        sstree = ss_tree.SS_Tree(seq, ss)
        for n in sstree:
            print n.data.sequence(), n.data.what()
            #print n.data.sequence() + "," + n.data.what() + "|",
        #for n in sstree:
        #    print n.data.ss

    def test_sub_ss_from_nodes(self):
        sstree = ss_tree.SS_Tree("GUG+CUC", "(.(+).)")
        nodes = sstree.tree.nodes[0:2]
        ss_data = sstree.sub_ss_from_nodes(nodes)
        #print ss_data.ss



def main():
    unittest.main()

if __name__ == '__main__':
    main()
