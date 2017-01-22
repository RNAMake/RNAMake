import unittest
from rnamake import motif, motif_tree, motif_state_tree, settings, exceptions
from rnamake import resource_manager
import rnamake.util as util
import build

class MotifStateTreeUnittest(unittest.TestCase):

    def setUp(self):
        self.rm = resource_manager.ResourceManager()

    def test_creation_from_mt(self):
        builder = build.BuildMotifTree(self.rm)
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(self.rm, mt)
        self.failUnless(len(mt) == len(mst), "did not build mst properly")

    def test_copy(self):
        mst = self._get_sub_tree()
        mst_copy = motif_state_tree.MotifStateTree.copy(mst)
        mst.add_connection(2, 3)

        self.failUnless(len(mst) == len(mst_copy))
        self.failUnless(mst.num_connections() == 1)
        self.failUnless(mst_copy.num_connections() == 0)

    def _get_sub_tree(self):
        mst = motif_state_tree.MotifStateTree(self.rm)
        mst.add_state(m_name="HELIX.IDEAL.2")
        mst.add_state(m_name="NWAY.1GID.0")
        mst.add_state(m_name="HELIX.IDEAL.2")
        mst.add_state(m_name="HELIX.IDEAL.2", parent_index=1)
        return mst

    def test_add_mst(self):
        mst = self._get_sub_tree()
        mst2 = self._get_sub_tree()

        mst.add_mst(mst2)
        self.failUnless(len(mst) == 8)

    def test_add_state(self):
        mst = motif_state_tree.MotifStateTree(self.rm)
        ms1 = self.rm.get_state(name="HELIX.IDEAL.2")
        ms2 = self.rm.get_state(name="HELIX.IDEAL.2")
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
            mst.add_state(ms2, parent_end_index=1, parent_end_name="A4-A5")

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
            mst.add_state(ms2, parent_end_name="A1-85")

        # invalid parent_end_name, cannot be found as an end in motif
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms2, parent_end_name="FAKE")

        mst.add_state(ms2)

    def test_add_connection(self):
        mst = motif_state_tree.MotifStateTree(self.rm)
        mst.add_state(m_name="HELIX.IDEAL.2")
        mst.add_state(m_name="NWAY.1GID.0")
        mst.add_state(m_name="HELIX.IDEAL.2")
        ms3 = self.rm.get_state(name="HELIX.IDEAL.2")

        # try connecting thru an already used end position
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_connection(1, 2, "A141-A162")

        mst.add_connection(1, 2)
        rna_struc = mst.get_structure()
        self.failUnless(len(rna_struc.chains()) == 1)
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_state(ms3, parent_end_index=1)

        self.failUnless(mst.add_state(ms3) == -1)
        with self.assertRaises(exceptions.MotifStateTreeException):
            mst.add_connection(1, 2)

    def test_remove_node(self):
        builder = build.BuildMotifTree(self.rm)
        mt = builder.build(3)
        mst = motif_state_tree.MotifStateTree(self.rm, mt)
        mst.remove_node(mst.last_node().index)
        self.failUnless(len(mst) != 2, "did not properly remove node")

    def test_remove_node_level(self):
        pass

    def test_replace_state(self):
        builder = build.BuildMotifTree(self.rm)
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(self.rm, mt)
        rstate = self.rm.get_state(name="HELIX.IDEAL.2")
        mst.replace_state(2, rstate)

        mt2 = motif_tree.MotifTree(self.rm)
        mt2.set_sterics(0)
        for n in mst:
            m = self.rm.get_motif(name=n.data.name,
                                  end_name=n.data.get_end(0).name)
            mt2.add_motif(m)

        i = len(mst)-1
        d1 = mst.get_node(i).data.get_end(1).d
        d2 = mt2.get_node(i).data.get_end(1).d
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