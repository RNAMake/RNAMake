import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree
import rnamake.motif_tree_state as motif_tree_state
import rnamake.motif_tree_precomputer as motif_tree_precomputer
import rnamake.settings as settings
import rnamake.x3dna as x3dna
import rnamake.motif_library as motif_library
import rnamake.motif_type as motif_type
import rnamake.motif_tree as motif_tree
import rnamake.util as util
import rnamake.motif as motif
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

    mtst = motif_tree_state.MotifTreeStateTree()
    mtst.add_state(h_mts)

    converter = motif_ensemble_tree.MTSTtoMETConverter()
    mtst2 = converter.convert(mtst, debug=1)
    mtst2.to_pdb('step.pdb')

def target_end_origin(hmotif):
    mt = motif_tree.MotifTree()
    mt.add_motif(hmotif)
    return mt.nodes[1].available_ends()[0].d()



me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
mts = me.motif_states[0].mts
mtst = motif_tree_state.MotifTreeStateTree()

motifs = []

helix_mlib = motif_library.MotifLibrary(motif_type.HELIX)
ideal = motif.str_to_motif(mts.build_string)
origin = target_end_origin(ideal)

mt = rnamake.motif_tree.MotifTree()
for i in range(1, 22):
    mtst.add_state(mts)
    m = mtst.to_pose()
    m.name = "HELIX.LE."+str(i)
    n = mt.add_motif(m, end_index=0, end_flip=0)
    m = mt.nodes[1].motif
    m.ends[0], m.ends[1] = m.ends[1], m.ends[0]
    #m.ends[0].flip()
    #m.ends[0].flipped = 0
    m.mtype = rnamake.motif_type.HELIX
    n_end_origin_1 = n.available_ends()[0].d()
    dist1 = util.distance(origin, n_end_origin_1)


    motifs.append(m)
    mt.remove_node_level()

#motifs[-1].to_pdb("test.pdb")
mtp = motif_tree_precomputer.MotifTreePrecomputer(name="HELIX",max_bps_per_end=0)
mtp.precompute_motifs(motifs)

me = motif_ensemble.MotifEnsemble("GC=GC", 0, 1)
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



