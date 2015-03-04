import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.motif_tree_state as motif_tree_state
import rnamake.settings as settings
import rnamake.x3dna as x3dna
import os
me = motif_ensemble.MotifEnsemble("GC=GC", 1, 0)
mts = me.motif_states[0].mts
print mts.name
name_elements = motif_tree_state.parse_db_name(mts.name)
print name_elements.start_index

mtst = motif_tree_state.MotifTreeStateTree()

path = settings.MOTIF_DIRS + "/helices/"
os.chdir(path)
x = x3dna.X3dna()
for i in range(22):
    mtst.add_state(mts)
    os.mkdir("HELIX.LE."+str(i))
    os.chdir("HELIX.LE."+str(i))
    m = mtst.to_pose()
    m.to_pdb("HELIX.LE."+str(i)+".pdb")
    x.get_basepairs("HELIX.LE."+str(i))
    os.chdir('..')




