import unittest
from rnamake import motif_connection, exceptions

class MotifConnectionUnittest(unittest.TestCase):

    def test_add_connection(self):
        mcs = motif_connection.MotifConnections()
        mcs.add_connection(1 , 2, "A1-A2", "B1-B2")
        self.failUnless(len(mcs) == 1)

    def test_remove_connection(self):
        mcs = motif_connection.MotifConnections()
        mcs.add_connection(1 , 2, "A1-A2", "B1-B2")

        with self.assertRaises(exceptions.MotifConnectionException):
            mcs.remove_connections_to(3)

        mcs.remove_connections_to(1)
        self.failUnless(len(mcs) == 0)

    def test_in_connection(self):
        mcs = motif_connection.MotifConnections()
        mcs.add_connection(1 , 2, "A1-A2", "B1-B2")

        self.failUnless(mcs.in_connection(1, "A1-A2") == 1)
        self.failUnless(mcs.in_connection(1, "A1-A3") == 0)
        self.failUnless(mcs.in_connection(2, "A1-A2") == 0)

def main():
    unittest.main()

if __name__ == '__main__':
    main()