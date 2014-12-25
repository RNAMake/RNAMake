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

