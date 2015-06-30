import unittest
import build
import rnamake.sqlite_library as sqlite_library
import rnamake.secondary_structure as secondary_structure
import rnamake.resource_manager as rm
import rnamake.motif_tree as motif_tree
import rnamake.motif_factory as motif_factory
import rnamake.ss_tree as ss_tree

class SecondaryStructureUnittest(unittest.TestCase):

    def test_assign_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        #for n in mt:
        #    print n.data.name

        p = mt.to_pose()
        #print p.secondary_structure()
        #print p.sequence()
        #mt.write_pdbs()

    def test_parse(self):
        m = rm.manager.get_motif("NWAY.1S72.18-01551-01634")
        ss = m.secondary_structure()
        seq = m.sequence()


        struct = secondary_structure.factory.get_structure(ss, seq)
        ss_chains = struct.reorient_ss_and_seq(0, 10)

        print struct.sequence()
        print struct.secondary_structure()



def main():
    unittest.main()

if __name__ == '__main__':
    main()