import unittest
import rnamake.ss_tree as ss_tree

class SSTreeUnittest(unittest.TestCase):

    def test_creation(self):
        #sstree = ss_tree.SS_Tree("GAG&CAG&CAC", "(.(&).(&).)")

        #sstree = ss_tree.SS_Tree("GUG+CUC", "(.(+).)")
        #exit()

        seq = "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA"
        ss  = "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...))))).."
        sstree = ss_tree.SS_Tree(seq, ss)
        for n in sstree:
            print n.data.sequence(), n.data.what()
            #print n.data.sequence() + "," + n.data.what() + "|",
        #for n in sstree:
        #    print n.data.ss

    def test_parse(self):
        seq = "UG&CA&CGACACAG"
        db  = "((&))&(......)"
        sstree = ss_tree.SS_Tree(seq, db)
        for n in sstree:
            print n.data.what(), n.data.sequence()



def main():
    unittest.main()

if __name__ == '__main__':
    main()
