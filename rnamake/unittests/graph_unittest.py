import unittest

from rnamake import graph, exceptions

class GraphUnittest(unittest.TestCase):

    #TODO check for weird uses of optional args in add_data
    def test_creation(self):
        g = graph.GraphDynamic()
        g.add_data(0)
        g.add_data(1)
        g.add_data(2, parent_index=0)
        g.add_data(3, parent_index=0)
        g.connect(2, 3)

        with self.assertRaises(exceptions.GraphIndexException):
            g.add_data(10, parent_index=10)


        g = graph.GraphStatic()
        # legend
        # N0 = Node 0
        # N0C0 = Connection 0 of Node 0
        # add first node with 2 possible connection
        g.add_data(0, n_children=2)
        #          N0
        #    N0C0 /  \  N0C1
        #        /    \
        #     None     None
        # add new node, connected to node 0. This new node: node 1 is connected
        # to node 0 using the connnection in the 0th position on both node 0 and
        # node 1. Do not need to specifiy parent_index=0 since adds to last node
        # without specifying
        # g.add_data(1, parent_index=0, parent_pos=0, child_pos=0, n_children=2)
        # means the same thing
        g.add_data(1, parent_pos=0, child_pos=0, n_children=2)
        #          N0
        #    N0C0 /  \  N0C1
        #    N1C0/    \
        #       N1   None
        #   N1C1|
        #       |
        #      None
        # add new node, connected to node 0.
        g.add_data(2, parent_index=0, parent_pos=1, child_pos=0, n_children=2)
        #          N0
        #    N0C0 /  \  N0C1
        #    N1C0/    \ N2C0
        #       N1    N2
        #   N1C1|      |N2C1
        #       |      |
        #      None   None
        g.connect(1, 2, 1, 1)
        #          N0
        #    N0C0 /  \  N0C1
        #    N1C0/    \ N2C0
        #       N1----N2
        #      N1C1 N2C1

        with self.assertRaises(exceptions.GraphInvalidEndException):
            g.connect(1, 2, 1, 1)

    def test_add_data(self):
        g = graph.GraphStatic()
        g.add_data(0, n_children=2)

        with self.assertRaises(exceptions.GraphIndexException):
            g.add_data(1, n_children=2, parent_index=10)

        with self.assertRaises(exceptions.GraphIndexException):
            g.add_data(1, n_children=2, index=0)

        with self.assertRaises(exceptions.GraphInvalidEndException):
            g.add_data(1, n_children=2, parent_pos=3)

        with self.assertRaises(exceptions.GraphInvalidEndException):
            g.add_data(1, n_children=2, parent_pos=-2)

        g.add_data(1, n_children=2, parent_pos=0, child_pos=1)

        with self.assertRaises(exceptions.GraphInvalidEndException):
            g.add_data(2, n_children=2, parent_index=0, parent_pos=0)

    def test_get_node(self):
        g = graph.GraphDynamic()
        g.add_data(0)
        g.add_data(1)

        self.failUnless(g.get_node(1).data == 1)

        with self.assertRaises(exceptions.GraphIndexException):
            g.get_node(10)

    def test_decrease_level(self):
        g = graph.GraphDynamic()

        #cant lower level past 0 which is default
        with self.assertRaises(exceptions.GraphException):
            g.decrease_level()

        g.increase_level()
        g.decrease_level()

        self.failUnless(g.level == 0)

    def test_removal(self):
        g = graph.GraphStatic()
        g.add_data(0, -1, -1, -1, 2)
        g.add_data(1, 0, 0, 0, 2)
        g.add_data(2, 1, 1, 0, 2)

        g.remove_node(1)
        g.add_data(3, 0, 0, 0, 2)
        g.connect(3, 2, 1, 1)

        #for n in graph.transverse_graph(g, 2):
        #    print n.data

        g = graph.GraphStatic()
        g.add_data(0, -1, -1, -1, 2)
        g.add_data(1, 0, 0, 0, 2)
        g.add_data(2, 1, 1, 0, 2)
        g.remove_node(1)

    def test_iter(self):
        g = graph.GraphStatic()
        for n in g:
            continue



def main():
    unittest.main()

if __name__ == '__main__':
    main()
