import unittest
import rnamake.prediction.setup
import rnamake.motif_type
import rnamake.motif_library
import rnamake.motif_tree_state
import util
import random

class SetupUnittest(unittest.TestCase):

    def test_aget_bp_setup(self):
        target = ["GC","GC"]
        mtype = rnamake.motif_type.HELIX
        mlib = rnamake.motif_library.MotifLibrary(mtype)
        mlib.load_all()

        rnamake.prediction.setup.get_bp_step_prediction_lib(target, mlib)

    def test_helix_motif_states(self):
        return
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath="test_helix.new.me")
        states = []
        for mts in mts_lib.motif_tree_states:
            name_elements = rnamake.motif_tree_state.parse_db_name(mts.name)
            if name_elements.flip_direction == 0 and name_elements.start_index == 0:
                states.append(mts)

        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        for mts in states:
            node = mtst.add_state(mts)
            if node is None:
                print mts.name
        mt = mtst.to_motiftree()

        for i in range(100):
            mtst = rnamake.motif_tree_state.MotifTreeStateTree()
            for j in range(100):
                mts = random.choice(states)
                node = mtst.add_state(mts)
                if node is None:
                    print mts.name
                    print mtst.last_node.mts.name
                    self.fail("failed")
            mt = mtst.to_motiftree()

    def test_error(self):
        return
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath="test_helix.new.me")
        mtst = rnamake.motif_tree_state.MotifTreeStateTree(sterics=0)
        mtst.add_state(mts_lib.get_state('helix.38-0-0-0-0-1-0'))
        node = mtst.add_state(mts_lib.get_state('helix.95-0-0-0-0-1-0'))
        mt = mtst.to_motiftree(sterics=0)



def main():
    unittest.main()

if __name__ == '__main__':
    main()
