import unittest
import rnamake.motif_tree_state
import rnamake.motif_type
import rnamake.settings
import rnamake.cluster
import util
import sys

class MotifTreeStateUnittest(unittest.TestCase):

    def test_creation(self):
        if util.UnittestState == util.UnittestType.BASIC:
            self.skipTest("test_creation is not a basic test")

        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        if len(mts_lib.motif_tree_states) != 510:
            self.fail("did not load all the mts")

        clusters = rnamake.cluster.cluster_mts(mts_lib.motif_tree_states,
                                               max_distance=0.1)
        if len(clusters) != 510:
            self.fail("no all mts are unique")

    def test_parse_db_name(self):
        name = "HELIX.LE.16-0-0-0-0-1-1"
        name_elements = rnamake.motif_tree_state.parse_db_name(name)
        if name_elements.motif_name != "HELIX.LE.16":
            self.fail()

def main():
    unittest.main()

if __name__ == '__main__':
    # TODO make this cleaner
    if len(sys.argv) > 1:
        try:
            util.UnittestState = int(sys.argv[1])
        except:
            pass
        sys.argv.pop()
    main()
