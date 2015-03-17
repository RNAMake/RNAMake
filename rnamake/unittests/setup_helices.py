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
import rnamake.basic_io as basic_io
import rnamake.motif_scorer as motif_scorer
import rnamake.util as util
import rnamake
import os


def motif_to_state(m, end_flip):
    mt = motif_tree.MotifTree()
    name = m.name+"-"+str(0)+"-"+str(0)+"-0-0-0-0"
    available_ends = [ m.ends[1] ]
    ends = [ None for e in m.ends]
    end_indexes = []
    for i, end in enumerate(m.ends):
        if 1 == i:
            continue
        ends[i] = end.state()

    beads = []
    m.get_beads([m.ends[0]])
    for b in m.beads:
        if b.btype != 0:
            beads.append(b.center)

    ms = motif_scorer.MotifScorer()
    score = ms.score(m)
    mts = motif_tree_state.MotifTreeState(name, 0, len(m.residues()), score, beads,
                         ends, end_flip, m.to_str())
    return mts




def write_mts_to_file(f, mts):
    f.write(mts.name + "|" + str(mts.score) + "|" + str(mts.size) + \
                "|" + str(mts.flip) +  "|" + mts.build_string + "|" + \
                basic_io.points_to_str(mts.beads) + "|")


    for i, end_state in enumerate(mts.end_states):
        if end_state is not None:
            f.write(end_state.to_str())
        else:
            f.write(" ")
        f.write("|")

    f.write("\n")
    f.flush()

def test_addition():
    me = motif_ensemble.MotifEnsemble("GC=GC", 0, 1)
    mts = me.motif_states[0].mts
    twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
    #helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
    helixs =  motif_tree_state.MotifTreeStateLibrary(libpath="HELIX.new.me")
    mtst = motif_tree_state.MotifTreeStateTree()
    mtst2 = motif_tree_state.MotifTreeStateTree()

    mtst.add_state(helixs.get_state('HELIX.LE.10-0-0-0-0-1-1'))

    for i in range(10):
        mtst2.add_state(mts)

    mtst.add_state(twoways.get_state('TWOWAY.2VQE.26-0-0-0-0-1-0'))
    node = mtst2.add_state(twoways.get_state('TWOWAY.2VQE.26-0-0-0-0-1-0'))
    print node
    mtst.to_pdb('test.pdb')
    mtst2.to_pdb('test2.pdb')
    mtst.nodes_to_pdbs()
    exit()



def check_answer(mts):

    print mts.name
    spl = mts.name.split("-")
    tail = "-".join(spl[1:])
    mts_lib = motif_tree_state.MotifTreeStateLibrary(libpath="HELIX.new.me")
    h_mts = mts_lib.get_state('HELIX.LE.10-'+tail)
    #h_mts = mts_lib.get_state('HELIX.LE.10-0-0-0-0-1-0')
    name_elements = motif_tree_state.parse_db_name(h_mts.name)
    m = rnamake.motif.str_to_motif(h_mts.build_string)
    m.to_pdb('LE.pdb')

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

f = open('HELIX_test.new.me','w')
for i in range(1, 22):
    break
    mtst.add_state(mts)
    m = mtst.to_pose()
    m.name = "HELIX.LE."+str(i)

    new_mts = motif_tree_state.motif_to_state(m, 0, 0)
    new_mts.name = m.name + "-0-0-0-0-1-0"

    last_node = mtst.last_node

    for j in range(len(new_mts.end_states)):
        if last_node.states[j] is None:
            continue
        new_mts.end_states[j] = last_node.states[j].copy()

    write_mts_to_file(f, new_mts)

me = motif_ensemble.MotifEnsemble("GC=GC", 0, 1)
mts = me.motif_states[0].mts
mtst = motif_tree_state.MotifTreeStateTree()
mt = motif_tree.MotifTree()


for i in range(1, 7):
    mtst.add_state(mts)
    m = mtst.to_pose()
    m.name = "HELIX.LE."+str(i)
    #mt.add_motif(m1,end_flip=0)
    #m = mt.nodes[1].motif

    dist = util.matrix_distance(mt.nodes[0].motif.ends[0].r(), m.ends[1].r())
    if dist < 0.1:
        m.ends[0], m.ends[1] = m.ends[1], m.ends[0]

    new_mts = motif_tree_state.motif_to_state(m, 0, 1)
    new_mts.flip = 1
    #new_mts.build_string = m.to_str()
    #new_mts.beads = beads
    new_mts.name = m.name + "-0-0-0-0-1-1"

    last_node = mtst.last_node
    m.to_pdb('h.'+str(i)+'.pdb')

    for j in range(len(new_mts.end_states)):
        if last_node.states[j] is None:
            continue
        print i
        print new_mts.end_states[j].r
        print last_node.states[j].r
        #new_mts.end_states[j] = last_node.states[j].copy()

    write_mts_to_file(f, new_mts)
    mt.remove_node_level()

f.close()

exit()

#motifs[-1].to_pdb("test.pdb")
mtp = motif_tree_precomputer.MotifTreePrecomputer(name="HELIX",max_bps_per_end=0)
mtp.precompute_motifs(motifs)

me = motif_ensemble.MotifEnsemble("GC=GC", 1, 0)
mts = me.motif_states[0].mts


test_addition()
#check_answer(mts)
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



