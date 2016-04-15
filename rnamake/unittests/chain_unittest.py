import unittest

from rnamake import chain, io, exceptions, settings

import util
import is_equal


class ChainUnittest(unittest.TestCase):

    def setUp(self):
        path = settings.UNITTEST_PATH + "resources/chain_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        chains = []
        for l in lines:
            chain = io.str_to_chain(l)
            chains.append(chain)

        self.chains = chains

    def test_creation(self):
        """creating a new object should never return an error"""

        try:
            c = chain.Chain()
        except:
            self.fail("was not expecting an error upon initation")

    def test_first_and_last(self):
        """test first and last functions, should return the first and last
           residues of the chain respectively.

           Also tests whether an empty chain will return exception when calling
           either
        """
        c = self.chains[0]
        if c.residues[0]  != c.first() or \
           c.residues[-1] != c.last():
            self.fail()

        chain_2 = chain.Chain()

        with self.assertRaises(exceptions.ChainException):
            chain_2.first()

        with self.assertRaises(exceptions.ChainException):
            chain_2.last()

    def test_subchain(self):
        c = self.chains[0]

        sub1 = c.subchain(0, 7)
        if len(sub1) != 7:
            self.fail()

        with self.assertRaises(exceptions.ChainException):
            c.subchain(-1, 7)

        with self.assertRaises(exceptions.ChainException):
            c.subchain(start_res=c.residues[10])

    def test_subchain_2(self):
        chain = self.chains[0]
        sub = chain.subchain(0)

        if not is_equal.are_chains_equal(chain, sub):
            self.fail("did not shallow copy chain correctly")

    def test_copy(self):
        c = self.chains[0]
        chain_copy = c.copy()
        if not is_equal.are_chains_equal(c, chain_copy):
            self.fail("did not copy chain correctly")

    def test_str(self):
        c = self.chains[0]
        s = c.to_str()
        c_copy = io.str_to_chain(s)

        if not is_equal.are_chains_equal(c, c_copy, check_uuid=0):
            self.fail("did not copy chain correctly")



def main():
    unittest.main()

if __name__ == '__main__':
    main()
