import unittest

from rnamake import graph, exceptions

class GraphUnittest(unittest.TestCase):

    def test_creation(self):
        g = graph.GraphDynamic()
        g.add_data(0)
        g.add_data(1)
        g.add_data(2, parent_index=0)
        g.add_data(3, parent_index=0)
        g.connect(2, 3)

        with self.assertRaises(exceptions.GraphException):
            g.add_data(10, parent_index=10)

        g = graph.GraphStatic()
        g.add_data(0, -1, -1, -1, 2)
        g.add_data(1, 0, 0, 0, 2)
        g.add_data(2, 0, 1, 0, 2)
        g.connect(1, 2, 1, 1)

    def test_get_node(self):
        g = graph.GraphDynamic()
        g.add_data(0)
        g.add_data(1)

        self.failUnless(g.get_node(1).data == 1)

        with self.assertRaises(exceptions.GraphException):
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
