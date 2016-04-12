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

        def get_chains(lines):
            chains = []
            for l in lines:
                chain = io.str_to_chain(l)
                chains.append(chain)
            return chains

        self.chains = util.supress_log_output(get_chains, lines)

    def test_creation(self):
        """
        creating a new object should never return an error
        """

        try:
            c = chain.Chain()
        except:
            self.fail("was not expecting an error upon initation")

    def test_first_and_last(self):
        c = self.chains[0]
        if c.residues[0]  != c.first() or \
           c.residues[-1] != c.last():
            self.fail()

        chain_2 = chain.Chain()

        with self.assertRaises(exceptions.ChainException):
            first = chain_2.first()

    def test_subchain(self):
        chain = self.chains[0]

        sub1 = chain.subchain(0, 7)
        if len(sub1) != 7:
            self.fail()

        try:
            sub2 = chain.subchain(-1, 7)
            self.fail()
        except ValueError:
            pass
        except:
            self.fail("unexpected error")

        sub3 = chain.subchain(0,-1)
        if len(sub3) != 156:
            self.fail()

    def test_subchain_2(self):
        chain = self.chains[0]
        sub = chain.subchain(0)

        if not is_equal.are_chains_equal(chain, sub):
            self.fail("did not shallow copy chain correctly")

    def test_copy(self):
        chain = self.chains[0]
        chain_copy = chain.copy()
        if not is_equal.are_chains_equal(chain, chain_copy):
            self.fail("did not copy chain correctly")






def main():
    unittest.main()

if __name__ == '__main__':
    main()
