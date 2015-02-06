import rnamake
import random

def test_add_state():
    mtype = rnamake.motif_type.TWOWAY
    mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(mtype)
    f = open("test_add_state.dat", "w")
    for j in range(100):
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        for i in range(100):
            mts = random.choice(mts_lib.motif_tree_states)
            mtst.add_state(mts)
        for n in mtst.nodes:
            f.write(n.mts.name + " ")
        f.write("\n")
    f.close()

test_add_state()


