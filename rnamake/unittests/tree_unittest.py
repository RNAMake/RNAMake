import unittest
import rnamake.tree as tree

class TreeUnittest(unittest.TestCase):

    def test_creation(self):
        t = tree.TreeDynamic()
        t.add_data(0)
        t.add_data(1)

        t.remove_node(index=1)

        for n in t:
            print n.index

        t = tree.TreeStatic()
        t.add_data(0, 2)
        t.add_data(1, 2, 0, 0)
        t.add_data(2, 2, 0, 1)





def main():
    unittest.main()

if __name__ == '__main__':
    main()
