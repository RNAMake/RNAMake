import rnamake.settings
import rnamake.resource_manager
import rnamake.motif_tree

rm = rnamake.resource_manager.ResourceManager()

def simple_mt():
    mt = rnamake.motif_tree.MotifTree()
    m = rm.get_motif("HELIX.IDEAL")
    mt.add_motif(m)
    mt.add_motif(m)
    return mt

def simple_mt_with_head():
    m = rm.get_motif("HELIX.IDEAL")
    mt = rnamake.motif_tree.MotifTree(m)
    mt.add_motif(m)
    mt.add_motif(m)
    return mt

def simple_mt_helix(size=10):
    m = rm.get_motif("HELIX.IDEAL")
    mt = rnamake.motif_tree.MotifTree()
    for i in range(0, size):
        mt.add_motif(m)
    return mt
