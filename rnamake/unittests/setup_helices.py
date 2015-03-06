import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.motif_tree_state as motif_tree_state
import rnamake.motif_tree_precomputer as motif_tree_precomputer
import rnamake.settings as settings
import rnamake.x3dna as x3dna
import rnamake
import os

def check_answer(mts):

    print mts.name
    spl = mts.name.split("-")
    tail = "-".join(spl[1:])
    mts_lib = motif_tree_state.MotifTreeStateLibrary(libpath="HELIX.new.me")
    h_mts = mts_lib.get_state('HELIX.LE.10-'+tail)
    name_elements = motif_tree_state.parse_db_name(h_mts.name)
    m = rnamake.motif.str_to_motif(h_mts.build_string)
    m.to_pdb('LE.pdb')

    m2 = rnamake.motif.str_to_motif(mts.build_string)
    m2.to_pdb("step.pdb")


me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
mts = me.motif_states[0].mts
mtst = motif_tree_state.MotifTreeStateTree()

motifs = []
mt = rnamake.motif_tree.MotifTree()
for i in range(22):
    mtst.add_state(mts)
    m = mtst.to_pose()
    m.name = "HELIX.LE."+str(i)
    mt.add_motif(m)
    #mt.nodes[1].motif.ends[0].flip()
    #mt.nodes[1].motif.ends[0].flipped=0
    m = mt.nodes[1].motif
    m.mtype = rnamake.motif_type.HELIX
    #m.ends[0], m.ends[1] = m.ends[1], m.ends[0]
    motifs.append(m)
    mt.remove_node_level()

#motifs[-1].to_pdb("test.pdb")
mtp = motif_tree_precomputer.MotifTreePrecomputer(name="HELIX",max_bps_per_end=0)

mtp.precompute_motifs(motifs)

me = motif_ensemble.MotifEnsemble("GC=GC", 1, 1)
mts = me.motif_states[0].mts


check_answer(mts)
exit()




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



