import unittest
import rnamake.ss_tree as ss_tree
import rnamake.motif_tree_topology as motif_tree_topology
import rnamake.resource_manager as resource_manager

class MotifTreeTopologyUnittest(unittest.TestCase):

    def test_creation(self):
        sstree = ss_tree.SS_Tree("(((+)))", "GAG+UUC")


        #sstree = ss_tree.SS_Tree("((.((+)).((+)).))", "GGAGG+CCAGG+CCACC")

        mtt = motif_tree_topology.MotifTreeTopology(sstree)
        for n in mtt:
            pass
            #for end in n.data.ends:
            #   print end.get_id(), resource_manager.manager.get_motif(end.get_id())

def main():
    unittest.main()

if __name__ == '__main__':
    main()
