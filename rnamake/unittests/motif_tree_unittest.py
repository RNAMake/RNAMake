import unittest
import build
import random
import rnamake.motif_tree as motif_tree
import rnamake.resource_manager as rm
import rnamake.secondary_structure_factory as ssfactory
import rnamake.motif_type as motif_type
from rnamake import exceptions, motif
import util

class MotifTreeUnittest(unittest.TestCase):

    def test_creation(self):
        mt = motif_tree.MotifTree()

    def test_add_motif_merger(self):
        mt = motif_tree.MotifTree()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m1)
        mt.add_motif(m2)
        if len(mt) != 2:
            self.fail("did not add motifs properly")

        if len(mt.merger.get_structure().residues()) != 14:
            self.fail("merger did not result in the right number of residues")

        if len(mt.merger.res_overrides.values()) != 2:
            self.fail("did not create the correct number of residue overrides")

        if len(mt.merger.bp_overrides.values()) != 1:
            self.fail("did not create the correct number of bp overrides")

    def test_add_motif_2(self):
        mt = motif_tree.MotifTree()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m1)

        #can never use parent_end_index=0 for a tree as that is where that node
        #is already connected to another node
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_end_index=0)

        #supplied parent_end_index and parent_end_name
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_end_index=1, parent_end_name="A1-A8")

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
            mt.add_motif(m2, parent_end_name="A4-A5")

        #invalid parent_end_name, cannot be found as an end in motif
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m2, parent_end_name="FAKE")

    def test_remove_node(self):
        mt = motif_tree.MotifTree()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m1)
        mt.add_motif(m2)
        mt.remove_node(1)
        if len(mt) != 1:
            self.fail("did not remove node correctly")

        if len(mt.merger.res_overrides.values()) != 0:
            self.fail(("did not remove residue overrides"))

        if len(mt.merger.bp_overrides.values()) != 0:
            self.fail(("did not remove bp overrides"))

    def test_remove_node_2(self):
        mt = motif_tree.MotifTree()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m3 = rm.manager.get_motif(name="HELIX.IDEAL.2")

        mt.add_motif(m1)
        mt.add_motif(m2)
        mt.add_motif(m3)

        mt.remove_node(1)
        if len(mt) != 2:
            self.fail("did not remove node correctly")

        if len(mt.merger.get_structure().chains()) != 4:
            print len(mt.merger.get_structure().chains())
            self.fail("did not get the correct number of chains")

        if len(mt.merger.get_structure().ends) != 4:
            self.fail("did not get the correct number of ends")

        #mt.merger.get_structure().to_pdb("test.pdb", renumber=1)

    def test_remove_node_3(self):
        mt = motif_tree.MotifTree()
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m)

        mt.remove_node(0)
        if len(mt) != 0:
            self.fail("did not remove node correctly")

    def test_remove_node_level(self):
        mt = motif_tree.MotifTree()
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

    def test_get_residues(self):
        mt = motif_tree.MotifTree()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m1)
        mt.add_motif(m2)
        res = []
        for n in mt:
            res.extend(n.data.residues())
        for r in res:
            if mt.merger.get_residue(r.uuid) is None:
                self.fail("could not find residue")

    def test_secondary_structure(self):
        mt = motif_tree.MotifTree()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m1)
        mt.add_motif(m2)
        ss = mt.secondary_structure()
        if ss.sequence() != "GGGGGGG&CCCCCCC":
            self.fail("did not get correct sequence")

        if ss.dot_bracket() != "(((((((&)))))))":
            self.fail("did not get correct dot bracket")

    def test_topology_to_str(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        s = mt.topology_to_str()
        mt2 = motif_tree.motif_tree_from_topology_str(s)
        if len(mt) != len(mt2):
            self.fail("did not get the right number of motifs")

        for i in range(len(mt)):
            n1 = mt.get_node(i)
            n2 = mt2.get_node(i)

            if n1.data.name != n2.data.name:
                self.fail("node " + str(i) + " did not have the same name")

            if n1.data.ends[0].name() != n2.data.ends[0].name():
                self.fail("node " + str(i) + " did not have the same end name")

    def test_replace_motif(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        m_new = builder.libs[0].get_random()
        #mt.to_pdb("test.pdb", renumber=1)
        #mt.write_pdbs("org")

        mt.replace_motif(2, m_new)
        #mt.write_pdbs()
        #mt.to_pdb("test_new.pdb", renumber=1)

    def test_pretty_print(self):
        mt = motif_tree.MotifTree();
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2");
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2");
        m3 = rm.manager.get_motif(name="HELIX.IDEAL.2");
        m4 = rm.manager.get_motif(name="HELIX.IDEAL.2");
        nway = rm.manager.get_motif(name="NWAY.1GID.0");

        mt.add_motif(m1);
        mt.add_motif(nway);
        mt.add_motif(m2);
        mt.add_motif(m3, 1);
        mt.add_motif(m4)


        #print mt.to_pretty_str()

    def _get_sub_motif_tree(self):
        mt = motif_tree.MotifTree();
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m3 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        nway = rm.manager.get_motif(name="NWAY.1GID.0")
        mt.add_motif(m1)
        mt.add_motif(nway)
        mt.add_motif(m2)
        mt.add_motif(m3, 1)

        return mt

    def test_print_pretty_2(self):
        mt = motif_tree.MotifTree();
        mt.option('sterics', 0)
        mt.add_motif_tree(self._get_sub_motif_tree())
        mt.add_motif_tree(self._get_sub_motif_tree())
        mt.add_motif_tree(self._get_sub_motif_tree(), parent_index=2,
                          parent_end_name="A1-A8")
        #print mt.to_pretty_str()
        #mt.write_pdbs()

    def test_get_build_points(self):
        mt = motif_tree.MotifTree()
        mt.add_motif(m_name="HELIX.IDEAL.2")

        build_points = mt.get_build_points()

        self.failUnless(len(build_points) == 1)
        self.failUnless(build_points[0].node.index == 0)
        self.failUnless(build_points[0].end_index == 1)

    def test_get_structure(self):
        mt = motif_tree.MotifTree()
        mt.add_motif(m_name="HELIX.IDEAL.2")
        mt.add_motif(m_name="HELIX.IDEAL.2")
        rna_struct = mt.get_structure()
        m = motif.Motif(r_struct=rna_struct)

        mt2 = motif_tree.MotifTree()
        mt2.add_motif(m)
        self.failUnless(len(mt2) == 1)

    def test_connections(self):
        mt = motif_tree.MotifTree()
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m3 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        nway = rm.manager.get_motif(name="NWAY.1GID.0")
        mt.add_motif(m1)
        mt.add_motif(nway)
        mt.add_motif(m2)

        # try connecting through 0th end position
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_connection(1, 2, "A138-A180")

        # try connecting thru an already used end position
        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_connection(1, 2, "A141-A162")

        mt.add_connection(1, 2)
        rna_struc = mt.get_structure()
        self.failUnless(len(rna_struc.chains()) == 1)

        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_motif(m3, parent_end_index=1)

        self.failUnless(mt.add_motif(m3) == -1)

        with self.assertRaises(exceptions.MotifTreeException):
            mt.add_connection(1, 2)

    def test_copy(self):
        mt = self._get_sub_motif_tree()
        mt_copy = mt.copy()
        mt.add_connection(2, 3)

        self.failUnless(len(mt) == len(mt_copy))
        self.failUnless(len(mt.connections.connections) == 1)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
