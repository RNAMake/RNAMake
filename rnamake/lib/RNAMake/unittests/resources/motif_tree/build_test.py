import rnamake.motif_library as motif_library
import rnamake.motif_tree as motif_tree
import random

def test_motif_tree_to_str():
    mlib = motif_library.unique_twoway_lib()
    f = open("test_motif_tree_to_str.dat", "w")
    for i in range(0, 100):
        mt = motif_tree.MotifTree()
        n = random.randint(10,20)
        j = 0
        count = 0
        while j < n:
            node = mt.add_motif(random.choice(mlib.motifs()), end_index=0)
            if node is not None:
                j += 1
            count += 1
            if count > 1000:
                break
        f.write(mt.to_str() + "\n")
    f.close()



test_motif_tree_to_str()
