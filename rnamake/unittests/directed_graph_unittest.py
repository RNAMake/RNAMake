import unittest

from rnamake import directed_graph, motif_directed_graph
from rnamake import resource_manager

class DirectedGraphUnittest(unittest.TestCase):

    def test_creation(self):
        dg = directed_graph.DirectedGraph()
        dg.add_node(99, 2)
        dg.add_node(97, 2, parent_index=0)
        dg.add_node(98, 2, parent_index=0, parent_edge_index=0)
        dg.add_node(96, 2)

        self.failUnless(dg.get_parent_index(0) is None)
        self.failUnless(dg.get_parent_index(1) == 0)
        self.failUnless(dg.get_parent_index(2) == 0)
        self.failUnless(dg.get_parent(1) == 99)

        self.failUnless(dg.get_roots() == [0, 3])

        self.failUnless(dg.are_nodes_connected(0, 1) == True)
        self.failUnless(dg.are_nodes_connected(0, 4) == False)

        dg.remove_node(0)
        self.failUnless(dg.get_parent_index(1) is None)

        dg.add_node(95, 2, parent_index=1)
        dg.remove_edge(1, 4)

    def test_iter(self):
        dg = directed_graph.DirectedGraph()
        dg.add_node(99, 3)
        dg.add_node(98, 2, parent_index=0, parent_edge_index=1)
        dg.add_node(97, 2)
        dg.add_node(96, 2, parent_index=1, parent_edge_index=1)
        dg.add_node(95, 2, parent_index=2, parent_edge_index=1)
        dg.add_edge(4, 0, 1, 0)

        nis = []
        for ni in dg:
            nis.append(ni)

        self.failUnless(nis == [0, 1, 3, 2, 4])

    def test_add_graph(self):
        dg = directed_graph.DirectedGraph()
        dg.add_node(99, 2)

        dg1 = directed_graph.DirectedGraph()
        dg1.add_node(100, 99)
        dg1.add_node(97, 99, parent_index=0)
        dg1.add_node(98, 99, parent_index=1)
        dg1.add_node(96, 99)
        dg1.add_edge(1, 3, 2, 2)
        dg1.add_edge(0, 3, 3, 3)

        dg.add_graph(dg1, 0, 0, 0, 0)
        self.failUnless(dg.get_parent_index(1) == 0)
        self.failUnless(dg.get_parent_index(4) is None)
        self.failUnless(dg.are_nodes_connected(1, 4))




class MotifDirectedGraphUnittest(unittest.TestCase):

    def test_creation(self):
        rm = resource_manager.ResourceManager()
        mdg = motif_directed_graph.MotifDirectedGraph()
        m = rm.get_motif(name="HELIX.IDEAL.2")

        print mdg.add_motif(m)




def main():
    unittest.main()

if __name__ == '__main__':
    main()
