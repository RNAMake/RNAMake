import unittest
import rnamake.graph as graph

class GraphUnittest(unittest.TestCase):

    def test_creation(self):
        g = graph.GraphDynamic()
        g.add_data(0)
        g.add_data(1)
        g.add_data(2, 0)
        g.add_data(3, 0)
        g.connect(2, 3)

        g = graph.GraphStatic()
        g.add_data(0, -1, -1, -1, 2)
        g.add_data(1, 0, 0, 0, 2)
        g.add_data(2, 0, 1, 0, 2)
        g.connect(1, 2, 1, 1)

    def test_removal(self):
        g = graph.GraphStatic()
        g.add_data(0, -1, -1, -1, 2)
        g.add_data(1, 0, 0, 0, 2)
        g.add_data(2, 1, 1, 0, 2)

        g.remove_node(1)
        g.add_data(3, 0, 0, 0, 2)
        g.connect(3, 2, 1, 1)

        for n in graph.transverse_graph(g, 2):
            print n.data

        g = graph.GraphStatic()
        g.add_data(0, -1, -1, -1, 2)
        g.add_data(1, 0, 0, 0, 2)
        g.add_data(2, 1, 1, 0, 2)
        g.remove_node(1)




def main():
    unittest.main()

if __name__ == '__main__':
    main()
