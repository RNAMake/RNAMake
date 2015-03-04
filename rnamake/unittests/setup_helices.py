import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.motif_tree_state as motif_tree_state
import rnamake.settings as settings
import rnamake.x3dna as x3dna
import os

me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
mts = me.motif_states[0].mts

mtst = motif_tree_state.MotifTreeStateTree()

path = settings.MOTIF_DIRS + "/helices/"
os.chdir(path)
x = x3dna.X3dna()
for i in range(22):
    mtst.add_state(mts)
    os.mkdir("HELIX.LE."+str(i))
    os.chdir("HELIX.LE."+str(i))
    mtst.to_pdb("HELIX.LE."+str(i)+".pdb")
    x.get_basepairs("HELIX.LE."+str(i))
    os.chdir('..')




