import rnamake.motif as motif
import rnamake.motif_tree as motif_tree
import rnamake.resource_manager as resource_manager

rm = resource_manager.ResourceManager()

ggaa_motif = motif.Motif("GGAA_tetraloop")
gaaa_motif = motif.Motif("GAAA_tetraloop")

mt = motif_tree.MotifTree()

mt.add_motif(rm.get_motif("AU=GC"))
mt.add_motif(ggaa_motif)

mt.write_pdbs()




