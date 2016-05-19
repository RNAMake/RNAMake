import unittest
from rnamake import motif_state_merger, motif_state_tree, motif_tree
from rnamake import resource_manager as rm
import build

class MotifStateMergerUnittest(unittest.TestCase):
    pass

    """
    def test_creation(self):
        ms = rm.manager.get_state(name="HELIX.IDEAL.2")
        mst = motif_state_tree.MotifStateTree()
        mst.add_state(ms)
        mst.add_state(ms)

        merger = motif_state_merger.MotifStateMerger()
        merger.merge(mst)
        rna_struct = merger._build_structure()
        ms_seq = rna_struct.sequence()

        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt =  motif_tree.MotifTree()
        mt.add_motif(m)
        mt.add_motif(m)

        m_seq = mt.secondary_structure().sequence()
        if ms_seq != m_seq:
            self.fail("did not produce correct sequence")

    def test_build_sequence(self):
        builder = build.BuildMotifTree()
        mt = builder.build(5)

        mst = motif_state_tree.MotifStateTree(mt)

        merger = motif_state_merger.MotifStateMerger()
        merger.merge(mst)
        rna_struct = merger._build_structure()
        ms_seq = rna_struct.sequence()
        m_seq = mt.secondary_structure().sequence()

        if ms_seq != m_seq:
            self.fail()

    def test_secondary_structure(self):
        ms = rm.manager.get_state(name="HELIX.IDEAL.2")
        mst = motif_state_tree.MotifStateTree()
        mst.add_state(ms)
        mst.add_state(ms)

        merger = motif_state_merger.MotifStateMerger()
        merger.merge(mst)
        ss = merger.secondary_structure()

    """

def main():
    unittest.main()

if __name__ == '__main__':
    main()