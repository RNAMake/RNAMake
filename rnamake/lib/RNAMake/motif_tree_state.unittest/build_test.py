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

def test_compare_last_node():
    mtype = rnamake.motif_type.TWOWAY
    mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(mtype)
    f = open("test_add_state.dat", "w")
    for j in range(100):
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        for i in range(10):
            mts = random.choice(mts_lib.motif_tree_states)
            mtst.add_state(mts)
        for n in mtst.nodes:
            f.write(n.mts.name + " ")
        avail_end = mtst.last_node.active_states()[0]
        f.write("|" + avail_ends.to_str() )
        f.write("\n")
    f.close()


#test_add_state()
test_compare_last_node()


