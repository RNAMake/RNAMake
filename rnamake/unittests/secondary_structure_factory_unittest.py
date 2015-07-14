import unittest
import rnamake.secondary_structure_factory as secondary_structure_factory
import rnamake.secondary_structure as secondary_structure
import rnamake.resource_manager as rm

class SecondaryStructureFactoryUnittest(unittest.TestCase):

    def test_creation(self):
        #m = rm.manager.get_motif("NWAY.1S72.18-01551-01634")
        seq = "GAGACAGACAC"
        db  = "(.(.).(.).)"

        ss = secondary_structure_factory.factory.get_structure(seq, db)
        mtt = ss.motif_topology_from_end(ss.ends[0])

    def test_from_motif(self):
        m = rm.manager.get_motif("HELIX.IDEAL.10")
        ss = secondary_structure.assign_secondary_structure(m)
        r = m.residues()[10]
        ss_r = ss.get_residue(uuid=r.uuid)
        if ss_r is None:
            raise ValueError("could not find corresponding ss res")





def main():
    unittest.main()

if __name__ == '__main__':
    main()