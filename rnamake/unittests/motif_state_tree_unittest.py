import unittest
from rnamake import motif, motif_tree, motif_state_tree, settings, exceptions
import rnamake.resource_manager as rm
import rnamake.util as util
import build

class MotifStateTreeUnittest(unittest.TestCase):


    def test_creation_from_mt(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(mt)
        if len(mt) != len(mst):
            self.fail("did not build mst properly")

    def test_copy(self):
        mst = self._get_sub_tree()
        mst_copy = mst.copy()
        mst.add_connection(2, 3)

        self.failUnless(len(mst) == len(mst_copy))
        self.failUnless(len(mst.connections.connections) == 1)


    def _get_sub_tree(self):
        mt = motif_tree.MotifTree()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m3 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        nway = rm.manager.get_motif(name="NWAY.1GID.0")
        mt.add_motif(m1)
        mt.add_motif(nway)
        mt.add_motif(m2)
        mt.add_motif(m3, 1)
        return mt

    def test_add_mst(self):
        mt = self._get_sub_tree()
        mt2 = self._get_sub_tree()

        mst = motif_state_tree.MotifStateTree(mt=mt)
        mst2 = motif_state_tree.MotifStateTree(mt=mt2)
        mst.add_mst(mst2)
        self.failUnless(len(mst) == 8)

    def test_add_state(self):
        mst = motif_state_tree.MotifStateTree()
        ms1 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms2 = rm.manager.get_state(name="HELIX.IDEAL.2")
        mst.add_state(ms1)

        # can only add the motif state with teh same unique indenitifer twice
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms1)

        # can never use parent_end_index=0 for a tree as that is where that node
        # is already connected to another node
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms2, parent_end_index=0)

        # supplied parent_end_index and parent_end_name
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms2, parent_end_index=1, parent_end_name="A1-A8")

        # must supply a motif or motif name
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state()

        # motif not found in resource manager
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(m_name="FAKE")

        # catches invalid parent_index
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms2, parent_index=2)

        # invalid parent_end_index, has only 0 and 1
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms2, parent_end_index=3)

        # invalid parent_end_name, is the name of end 0
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms2, parent_end_name="A4-A5")

        # invalid parent_end_name, cannot be found as an end in motif
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms2, parent_end_name="FAKE")

    def test_add_connection(self):
        mst = motif_state_tree.MotifStateTree()
        ms1 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms2 = rm.manager.get_state(name="HELIX.IDEAL.2")
        ms3 = rm.manager.get_state(name="HELIX.IDEAL.2")
        nway = rm.manager.get_state(name="NWAY.1GID.0")
        mst.add_state(ms1)
        mst.add_state(nway)
        mst.add_state(ms2)

        # try connecting through 0th end position
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_connection(1, 2, "A138-A180")

        # try connecting thru an already used end position
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_connection(1, 2, "A138-A180")

        mst.add_connection(1, 2)
        rna_struc = mst.get_structure()
        self.failUnless(len(rna_struc.chains()) == 1)
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms3, parent_end_index=1)

        self.failUnless(mst.add_state(ms3) == -1)
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_connection(1, 2)

    def _test_topology_to_str(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        mst = motif_state_tree.MotifStateTree(mt)
        s = mst.topology_to_str()
        mst2 = motif_state_tree.str_to_motif_state_tree(s)
        if len(mst2) != len(mt):
            self.fail("did not reconstituate motif_state_tree from topology")

    def test_remove_node(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mst = motif_state_tree.MotifStateTree(mt)
        mst.remove_node(mst.last_node().index)
        if len(mst) != 2:
            self.fail("did not properly remove node")

    def test_remove_node_level(self):
        pass

    def test_replace_state(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(mt)
        rstate = rm.manager.get_state(name="HELIX.IDEAL.2")
        mst.replace_state(2, rstate)

        mt2 = motif_tree.MotifTree()
        mt2.option('sterics', 0)
        for n in mst:
            m = rm.manager.get_motif(name=n.data.cur_state.name,
                                     end_name=n.data.cur_state.end_names[0])
            mt2.add_motif(m)

        i = len(mst)-1
        d1 = mst.get_node(i).data.cur_state.end_states[1].d
        d2 = mt2.get_node(i).data.ends[1].d()
        if util.distance(d1, d2) > 0.1:
            self.fail("did not properly replace state")

    def test_pretty_str(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mst = motif_state_tree.MotifStateTree(mt)
        s = mst.to_pretty_str()


def main():
    unittest.main()

if __name__ == '__main__':
    main()