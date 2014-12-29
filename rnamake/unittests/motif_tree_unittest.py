import unittest
import rnamake.motif_tree
import rnamake.resource_manager
import util
import instance

class MotifTreeUnittest(unittest.TestCase):

    def test_creation(self):
        mt = rnamake.motif_tree.MotifTree()
        sterics = mt.option('sterics')
        if sterics != 1:
            self.fail("cannot retreive correct value")

        mt = rnamake.motif_tree.MotifTree(sterics=0)
        sterics = mt.option('sterics')
        if sterics != 0:
            self.fail("cannot retreive correct value")

    def test_add_motif(self):
        mt = rnamake.motif_tree.MotifTree()
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL")
        mt.add_motif(m)
        mt.add_motif(m)
        if len(mt.nodes) != 3:
            self.fail("did not add motifs properly")

        node = mt.add_motif(m, end_index=1, end_flip=0)
        if node is not None:
            self.fail("did not perform sterics correctly")

    def test_remove_node(self):
        mt = rnamake.motif_tree.MotifTree()
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL")
        mt.add_motif(m)
        mt.add_motif(m)

        mt.remove_node(mt.last_node)
        if len(mt.nodes) != 2:
            self.fail("did not remove node correctly")

    def test_to_str_node(self):
        mt = instance.simple_mt()
        s = mt.to_str()
        mt2 = rnamake.motif_tree.str_to_motif_tree(s)



def main():
    unittest.main()

if __name__ == '__main__':
    main()
