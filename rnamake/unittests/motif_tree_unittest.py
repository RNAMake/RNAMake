import unittest
import build
import random
import rnamake.motif_tree as motif_tree
import rnamake.resource_manager as rm
import rnamake.secondary_structure_factory as ssfactory
import rnamake.motif_type as motif_type
import rnamake.eternabot.sequence_designer as sequence_designer
import rnamake.motif_tree_topology as motif_tree_topology
import util

class MotifTreeUnittest(unittest.TestCase):

    def test_creation(self):
        mt = motif_tree.MotifTree()

    def test_add_motif(self):
        mt = motif_tree.MotifTree()
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)
        if len(mt) != 2:
            self.fail("did not add motifs properly")

        if len(mt.merger.get_structure().residues()) != 14:
            self.fail("merger did not result in the right number of residues")

    def test_remove_node(self):
        mt = motif_tree.MotifTree()
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)

        mt.remove_node(1)
        if len(mt) != 1:
            self.fail("did not remove node correctly")

    def test_remove_node_2(self):
        mt = motif_tree.MotifTree()
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)
        mt.add_motif(m)

        mt.remove_node(1)
        if len(mt) != 2:
            self.fail("did not remove node correctly")

        if len(mt.merger.get_structure().chains()) != 4:
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

        #mt.merger.get_structure().to_pdb("test.pdb", renumber=1)

    def test_secondary_structure(self):
        rm.manager.add_motif("resources/motifs/tetraloop_receptor_min")
        mt = motif_tree.MotifTree()
        mt.add_motif(rm.manager.get_motif(name="tetraloop_receptor_min",
                                          end_name="A228-A246"))
        mt.add_motif(rm.manager.get_motif(name="HELIX.IDEAL.20"), parent_end_name="A221-A252")
        #mt.option('sterics', 0)
        mt.add_motif(rm.manager.get_motif(name="tetraloop_receptor_min",
                                          end_name="A221-A252"))

        ss = mt.secondary_structure()






def main():
    unittest.main()

if __name__ == '__main__':
    main()
