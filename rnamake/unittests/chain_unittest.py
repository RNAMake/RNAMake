import unittest
import rnamake.chain
import rnamake.io
import util


class ChainUnittest(unittest.TestCase):

    def setUp(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/chain_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        def get_chains(lines):
            chains = []
            for l in lines:
                chain = rnamake.io.str_to_chain(l)
                chains.append(chain)
            return chains

        self.chains = util.supress_log_output(get_chains, lines)

    def test_creation(self):
        try:
            c = rnamake.chain.Chain()
        except:
            self.fail("was not expecting an error upon initation")

    def test_str_to_chain(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/chain_strs.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        def get_chains(lines):
            chains = []
            for l in lines:
                chain = rnamake.io.str_to_chain(l)
                chains.append(chain)
            return chains

        chains = util.supress_log_output(get_chains, lines)
        if len(chains) == 0:
            self.fail("did not load chains")

    def test_first_and_last(self):
        chain = self.chains[0]
        if chain.residues[0]  != chain.first() or \
           chain.residues[-1] != chain.last():
            self.fail()

    def test_repr(self):
        try:
            chain = rnamake.chain.Chain()
        except:
            self.fail("chain __repr__ yeilded unexpected error")

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




def main():
    unittest.main()

if __name__ == '__main__':
    main()
