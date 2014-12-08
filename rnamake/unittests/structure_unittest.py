import unittest
import rnamake.structure
import rnamake.settings
import util

class StructureUnittest(unittest.TestCase):

    def test_build_chains(self):
        path = rnamake.settings.UNITTEST_PATH + "resources/p4p6.pdb"
        struct = util.supress_log_output(rnamake.structure.Structure,
                                            path)

        print "chains ", len(struct.chains)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
