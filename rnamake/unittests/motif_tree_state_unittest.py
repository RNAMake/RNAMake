import unittest
import rnamake.motif_tree_state
import rnamake.motif_type
import rnamake.settings
import rnamake.cluster
import rnamake.basepair
import util
import numerical
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

    def test_node_creation(self):
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        node = rnamake.motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])

    def test_node_copy(self):
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        node = rnamake.motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])
        cnode = node.copy()
        cnode.state.d += [0,0,-5]
        if numerical.are_points_equal(cnode.state.d, node.state.d):
            self.fail("copy was not sucessful")

    def test_aligner(self):
        return
        path = rnamake.settings.UNITTEST_PATH + \
               "/resources/test_node_align.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        node1 = rnamake.motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])
        node2 = rnamake.motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])

        aligner = rnamake.motif_tree_state.MotifTreeStateNodeAligner()
        for l in lines:
            spl = l.split(",")
            bp1 = rnamake.basepair.str_to_basepairstate(spl[0])
            bp2 = rnamake.basepair.str_to_basepairstate(spl[1])
            bp3 = rnamake.basepair.str_to_basepairstate(spl[2])
            node1.state = bp1
            node2.state = bp2
            aligner.transform_state(node1, node2)
            print node2.state.r
            print bp3.r
            exit ()

    def test_tree(self):
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()

        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        mtst = rnamake.motif_tree_state.MotifTreeStateTree(mts)

    def test_tree_add(self):
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]

        node = mtst.add_state(mts)
        if node is None or len(mtst.nodes) != 2:
            self.fail("did not add node properly")

    def test_tree_to_motiftree(self):
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        mt = mtst.to_motiftree()
        if mt.nodes[0].motif.name != "start":
            self.fail("did not convert properly")


        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        mtst = rnamake.motif_tree_state.MotifTreeStateTree(mts)
        mt = mtst.to_motiftree()

        for mts in mts_lib.motif_tree_states:
            node = mtst.add_state(mts)
            if node is None:
                continue
            mt = mtst.to_motiftree()
            return



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
