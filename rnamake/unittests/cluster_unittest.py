import unittest
import rnamake.cluster
import rnamake.motif_tree_state
import copy

class ClusterUnittest(unittest.TestCase):

    def test_cluster_mts(self):
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        mtss = [mts, copy.deepcopy(mts), copy.deepcopy(mts)]
        clusters = rnamake.cluster.cluster_mts(mtss)
        if len(clusters) != 1:
            self.fail("did not cluster properly")

        print mtss[1].end_state.d
        mtss[1].end_state.d -= [0,0,5]
        print mtss[1].end_state.d
        clusters = rnamake.cluster.cluster_mts(mtss)
        if len(clusters) != 2:
            self.fail("did not cluster properly")



def main():
    unittest.main()

if __name__ == '__main__':
    main()
