import unittest
import rnamake.settings as settings
import rnamake.ss_tree as ss_tree
import rnamake.secondary_structure as secondary_structure
import rnamake.motif_factory as motif_factory
import rnamake.secondary_structure_factory as sf

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

    def test_parse_tc(self):
        path = settings.MOTIF_DIRS + "/tertiary_contacts/TC.1DUQ.0"
        tc = motif_factory.factory.motif_from_file(path)
        #ss = secondary_structure.assign_secondary_structure(tc)
        #motif_factory.factory._setup_secondary_structure(tc)
        print tc.secondary_structure
        ss = sf.factory.get_structure(base_ss=tc.secondary_structure)
        print ss.ends


def main():
    unittest.main()

if __name__ == '__main__':
    main()
