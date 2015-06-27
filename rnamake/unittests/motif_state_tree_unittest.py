import unittest
import rnamake.motif_tree as motif_tree
import rnamake.motif_factory as motif_factory
import rnamake.motif_state_tree as motif_state_tree
import rnamake.resource_manager as rm
import rnamake.util as util
import build

class MotifStateTreeUnittest(unittest.TestCase):

    def _test_creation(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)

        mst = motif_state_tree.MotifStateTree()

        for i, n in enumerate(mt):
            if i == 0:
                mst.add_state(rm.manager.get_state(n.data.name))
            else:
                parent_index = 1000
                parent_end_index = -1
                for c in n.connections:
                    if c is None:
                        continue
                    if c.partner(n.index).index < parent_index:
                        parent_index =  c.partner(n.index).index
                        parent_end_index = c.end_index(c.partner(n.index).index)

                if parent_index == 1000:
                    raise ValueError("did not convert motif tree to motif state tree properly")
                mst.add_state(rm.manager.get_state(n.data.name), parent_index, parent_end_index)

    def _test_creation_from_mt(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(mt)

    def _test_to_mt(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mt.write_pdbs()

        for n in mt:
            print n.data.name

        mst = motif_state_tree.MotifStateTree(mt)

        mt2 = mst.to_motif_tree()
        mt2.to_pdb("test2.pdb")

    def test_replace_state(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(mt)
        rstate = rm.manager.get_state("HELIX.IDEAL.2")
        mst.replace_state(4, rstate)

        mt2 = motif_tree.MotifTree()
        for n in mst:
            m = rm.manager.get_motif(n.data.ref_state.name)
            mt2.add_motif(m)

        i = len(mst)-1
        d1 = mst.get_node(i).data.cur_state.end_states[1].d
        d2 = mt2.get_node(i).data.ends[1].d()

        if util.distance(d1, d2) > 0.1:
            self.fail("did not properly replace state")


    def test_specific(self):
        str = """HELIX.IDEAL.11
TWOWAY.2VQE.7-A416-A427
HELIX.IDEAL.16
TWOWAY.2VQE.13-A662-A743
HELIX.IDEAL.19
TWOWAY.1S72.29-01312-01342
HELIX.IDEAL.13
TWOWAY.1D4R.1-A17-B12
HELIX.IDEAL.3
TWOWAY.2VQE.9-A455-A477"""
        motif_names = str.split("\n")
        mt = motif_tree.MotifTree()
        for m_n in motif_names:
            mt.add_motif(rm.manager.get_motif(m_n))

        mst = motif_state_tree.MotifStateTree(mt)
        mt2 = mst.to_motif_tree()







def main():
    unittest.main()

if __name__ == '__main__':
    main()