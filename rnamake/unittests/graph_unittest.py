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




def main():
    unittest.main()

if __name__ == '__main__':
    main()
