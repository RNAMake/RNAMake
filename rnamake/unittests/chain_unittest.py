import unittest

from rnamake.chain import Chain
from rnamake import exceptions, settings, residue_type, motif_state

import util
import is_equal


class ChainUnittest(unittest.TestCase):

    def setUp(self):
        self.rts = residue_type.ResidueTypeSet()

        path = settings.UNITTEST_PATH + "resources/chain_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        chains = []
        for l in lines:
            c = Chain.from_str(l, self.rts)
            chains.append(c)

        self.chains = chains

    def test_creation(self):
        """creating a new object should never return an error"""

        try:
            c = Chain()
        except:
            self.fail("was not expecting an error upon initation")

    def test_first_and_last(self):
        """test first and last functions, should return the first and last
           residues of the chain respectively.

           Also tests whether an empty chain will return exception when calling
           either
        """
        c = self.chains[0]
        if c.residue(0)  != c.first() or \
           c.residue(-1) != c.last():
            self.fail()

        chain_2 = Chain()

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
            c.subchain(start_res=c.residue(10))

    def test_subchain_2(self):
        chain = self.chains[0]
        sub = chain.subchain(0)

        self.failUnless(is_equal.are_chains_equal(chain, sub))

    def test_copy(self):
        c = self.chains[0]
        c_copy = Chain.copy(c)
        if not is_equal.are_chains_equal(c, c_copy):
            self.fail("did not copy chain correctly")

    def test_str(self):
        c = self.chains[0]
        s = c.to_str()
        c_copy = Chain.copy(c)

        if not is_equal.are_chains_equal(c, c_copy, check_uuid=0):
            self.fail("did not copy chain correctly")

    def test_get_state(self):
        c = self.chains[0]
        cs = c.get_state()

        self.failUnless(len(c) == len(cs))

        cs_copy = motif_state.Chain.copy(cs)

        self.failUnless(len(cs) == len(cs_copy))
        self.failUnless(cs.first() == cs_copy.first())
        self.failUnless(cs.last() == cs_copy.last())

        s = cs.to_str()
        cs_copy = motif_state.Chain.from_str(s)

        self.failUnless(len(cs) == len(cs_copy))

def main():
    unittest.main()

if __name__ == '__main__':
    main()
2