import unittest
import rnamake.secondary_structure_factory as secondary_structure_factory
import rnamake.resource_manager as rm

class SecondaryStructureFactoryUnittest(unittest.TestCase):

    def test_creation(self):
        #m = rm.manager.get_motif("NWAY.1S72.18-01551-01634")
        seq = "GAGACAGACAC"
        db  = "(.(.).(.).)"

        ss = secondary_structure_factory.factory.get_structure(seq, db)
        print ss.ends[0].res1.num, ss.ends[0].res2.num
        mtt = ss.motif_topology_from_end(ss.ends[0])
        for n in mtt:
            print n



def main():
    unittest.main()

if __name__ == '__main__':
    main()