import unittest
import build
import rnamake.motif_tree as motif_tree
import rnamake.motif_type as motif_type
from rnamake import exceptions, motif, resource_manager

class MotifTreeUnittest(unittest.TestCase):

    def setUp(self):
        self.rm = resource_manager.ResourceManager()

    def test_creation(self):
        motif_tree.MotifTree(self.rm)

    def test_add_motif(self):
        mt = motif_tree.MotifTree(self.rm)
        m1 = self.rm.get_motif(name="HELIX.IDEAL.2")
        m2 = self.rm.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m1)

        self.failUnless(len(mt) == 1)

        #can never use parent_end_index=0 for a tree as that is where that node
        #is already connected to another node
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_end_index=0)

        #supplied parent_end_index and parent_end_name
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_end_index=1, parent_end_name="A4-A5")

        #must supply a motif or motif name
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif()

        #motif not found in resource manager
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m_name="FAKE")

        #catches invalid parent_index
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_index=2)

        #invalid parent_end_index, has only 0 and 1
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_end_index=3)

        #invalid parent_end_name, is the name of end 0
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_end_name="A1-A8")

        #invalid parent_end_name, cannot be found as an end in motif
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_end_name="FAKE")

    def test_remove_node(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.remove_node(1)
        if len(mt) != 1:
            self.fail("did not remove node correctly")

    def test_remove_node_2(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.remove_node(1)
        if len(mt) != 2:
            self.fail("did not remove node correctly")

        rna_struc = mt.get_structure()

        if rna_struc.num_chains() != 4:
            self.fail("did not get the correct number of chains")

        if rna_struc.num_ends()  != 4:
            self.fail("did not get the correct number of ends")

    def test_remove_node_3(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.remove_node(0)
        if len(mt) != 0:
            self.fail("did not remove node correctly")

    def test_remove_node_level(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.increase_level()
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2")

        mt.remove_node_level()
        self.failUnless(len(mt) == 1)

        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2")

        mt.remove_node_level(level=0)
        self.failUnless(len(mt) == 0)

    def test_secondary_structure(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2")

        ss = mt.secondary_structure()
        if ss.sequence() != "CCCCCCC&GGGGGGG":
            self.fail("did not get correct sequence")

        if ss.dot_bracket() != "(((((((&)))))))":
            self.fail("did not get correct dot bracket")

    def test_replace_motif(self):
        builder = build.BuildMotifTree(self.rm)
        mt = builder.build(10)
        m_new = builder.libs[0].get_random()

        mt.replace_motif(2, m_new)

        mt2 = motif_tree.MotifTree(self.rm)

    def _get_sub_motif_tree(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="NWAY.1GID.0")
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2", parent_index=1)

        return mt

    def test_add_motif_tree(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.set_sterics(0)
        mt.add_motif_tree(self._get_sub_motif_tree())
        mt.add_motif_tree(self._get_sub_motif_tree())
        mt.add_motif_tree(self._get_sub_motif_tree(), parent_index=2,
                          parent_end_name="A4-A5")
        #print mt.to_pretty_str()
        #mt.write_pdbs()

    def test_get_build_points(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")

        build_points = mt.get_build_points()

        self.failUnless(len(build_points) == 1)
        self.failUnless(build_points[0].node.index == 0)
        self.failUnless(build_points[0].end_index == 1)

    def test_get_structure(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2")
        rna_struct = mt.get_structure()

        self.rm.add_motif(rna_struct, motif_type.HELIX, "double_helix")
        m = self.rm.get_motif(name="double_helix")

        mt2 = motif_tree.MotifTree(self.rm)
        mt2.add_motif(m)
        self.failUnless(len(mt2) == 1)

    def test_connections(self):

        m3 = self.rm.get_motif(name="HELIX.IDEAL.2")
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="NWAY.1GID.0")
        mt.add_motif(m_name="HELIX.IDEAL.2")

        # try connecting thru an already used end position
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_connection(1, 2, "A141-A162")

        mt.add_connection(1, 2)
        rna_struc = mt.get_structure()
        self.failUnless(rna_struc.num_chains() == 1)

        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m3, parent_end_index=1)

        self.failUnless(mt.add_motif(m3) == -1)

        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_connection(1, 2)

    def test_copy(self):
        mt = self._get_sub_motif_tree()
        mt_copy = motif_tree.MotifTree.copy(mt)
        mt.add_connection(2, 3)

        self.failUnless(len(mt) == len(mt_copy))
        self.failUnless(len(mt._connections.connections) == 1)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
