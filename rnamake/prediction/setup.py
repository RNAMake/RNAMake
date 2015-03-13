import rnamake.motif as motif
import rnamake.motif_tree as motif_tree
import rnamake.util as util
import rnamake.cluster as cluster
import rnamake.motif_tree_precomputer as motif_tree_precomputer
import rnamake.settings as settings
import rnamake.motif_library as motif_library
import rnamake.motif_type as motif_type
import rnamake.motif_outputer as motif_outputer
import copy
import math

def bp_match_target(bp, target):
    bpstr = bp.res1.rtype.name[0]+bp.res2.rtype.name[0]
    return bpstr == target


def residue_in_matched_bp(m, r, target):
    bps = m.get_basepair(uuid1=r.uuid)
    for bp in bps:
        if r != bp.res1:
            bp.res1, bp.res2 = bp.res2, bp.res1
        if bp.bp_type == "cW-W" and bp_match_target(bp, target):
            return bp


def motif_from_bps(bps):
    # bps = [copy.deepcopy(bp) for bp in bps]
    m = motif.Motif()
    res = []
    for bp in bps:
        res.extend(bp.residues())
    m.structure._build_chains(res)
    m.basepairs = bps
    m.structure._cache_coords()
    m._cache_basepair_frames()
    m.ends = [bps[0], bps[-1]]
    return m


def target_end_origin(hmotif):
    mt = motif_tree.MotifTree()
    mt.add_motif(hmotif)
    return mt.nodes[1].available_ends()[0].d()

def get_bp_step_prediction_lib(targets, helix_mlib):
    not_seen = {}
    aligned_motifs = []
    mt = motif_tree.MotifTree()
    ideal = helix_mlib.get_motif("HELIX.IDEAL")
    origin = target_end_origin(ideal)
    name = "=".join(targets)
    base_dir = settings.RESOURCES_PATH + "prediction/"

    for m in helix_mlib.motifs():
        spl = m.name.split(".")
        #if spl[1] not in include:
        #    continue
        if spl[1] == "IDEAL" or spl[1] == "LE":
            continue
        for c in m.chains():
            for i in range(len(c.residues)-1):
                bp1 = residue_in_matched_bp(m, c.residues[i], targets[0])
                bp2 = residue_in_matched_bp(m, c.residues[i+1], targets[1])
                if bp1 is None or bp2 is None:
                    continue
                m = motif_from_bps((bp1, bp2))
                n = mt.add_motif(m)
                n_end_origin = n.available_ends()[0].d()
                if util.distance(n_end_origin, origin) > 3:
                    mt.remove_node(mt.last_node)
                    mt.add_motif(m, end_flip=1)
                    mt.nodes[1].motif.ends[0].flipped = 0
                aligned_motifs.append(mt.nodes[1].motif)
                mt.remove_node(mt.last_node)

    kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
    kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)
    clusters = cluster.cluster_motifs(aligned_motifs)
    motifs = []
    nmotifs = len(aligned_motifs)
    test_helix = helix_mlib.get_motif("HELIX.IDEAL.2")
    f = open(base_dir + "/" + name + ".pop", "w")
    for i, c in enumerate(clusters):
        mt.add_motif(c.motifs[0])
        node = mt.add_motif(test_helix, end_index=0, end_flip=0)
        if mt.nodes[1].motif.ends[1].r()[1][0] > 0:
            mt.nodes[1].motif.ends[1].flip()
            mt.nodes[1].motif.ends[1].flipped=0
        pop = float(len(c.motifs)) / float(nmotifs)
        energy = -kBT*math.log(pop)
        mt.nodes[1].motif.name = name+"."+str(i)
        f.write(name+"."+str(i) + " " + str(len(c.motifs)) + " " + \
                str(pop) + " " + str(energy) + "\n")
        motifs.append(mt.nodes[1].motif)
        mt.remove_node_level()
    f.close()
    #mo = motif_outputer.MotifOutputer()
    #for m in motifs:
    #    mo.add_motif(m,0)
    #mo.to_pdb()
    #exit()
    #for i, m in enumerate(motifs):
    #    mt.add_motif(m, end_index=1)
    #    mt.nodes[1].motif.to_pdb("motif."+str(i)+".pdb")
    #    mt.remove_node(mt.last_node)
    mtp = motif_tree_precomputer.MotifTreePrecomputer(name=base_dir+"/"+name,
            max_bps_per_end=0)
    mtp.precompute_motifs(motifs)


class SSandSeqCluster(object):
    def __init__(self, motif):
        self.ss = motif.secondary_structure()
        self.seq = motif.sequence()
        self.motifs = [ motif ]

    def same_seq_and_ss(self, m):
        if len(self.ss) != len(m.secondary_structure()):
            return 0

        if self.ss == m.secondary_structure() and self.seq == m.sequence():
            return 1

        ss = m.secondary_structure()
        seq = m.sequence()
        ss_spl = ss.split('&')
        seq_spl = seq.split('&')

        rseq = seq_spl[1]+'&'+seq_spl[0]
        new_ss = ss_spl[1]+'&'+ss_spl[0]
        rss = ""
        for e in new_ss:
            if   e == '(':
                rss += ')'
            elif e == ')':
                rss += '('
            else:
                rss += e

        if self.ss == rss and self.seq == rseq:
            return 1

        else:
            return 0


def get_ss_name(ss):
    name = ""
    for e in ss:
        if e == '(' or e == ')':
            name += 'P'
        elif e == '.':
            name += 'U'
        elif e == '&':
            name += '-'
    return name

def get_two_way_prediction_lib(c, origin, test_helix):
    mt = motif_tree.MotifTree()

    i = 0
    aligned_motifs = []
    for m in c.motifs:
        n = mt.add_motif(m, end_flip=0)
        n_end_origin_1 = n.available_ends()[0].d()
        mt.remove_node_level()
        n = mt.add_motif(m, end_flip=1)
        n_end_origin_2 = n.available_ends()[0].d()

        dist1 = util.distance(origin, n_end_origin_1)
        dist2 = util.distance(origin, n_end_origin_2)
        if dist2 > dist1:
            aligned_motifs.append(mt.nodes[1].motif)
        else:
            mt.remove_node_level()
            mt.add_motif(m, end_flip=0)
            mt.nodes[1].motif.ends[0].flipped = 0
            aligned_motifs.append(mt.nodes[1].motif)

        mt.remove_node_level()
        i += 1

    kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
    kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)
    clusters = cluster.cluster_motifs(aligned_motifs)
    nmotifs = len(aligned_motifs)
    motifs = []

    base_dir = settings.RESOURCES_PATH + "prediction/pdb_ensembles/"
    seq_spl = c.seq.split('&')
    seq_name = '-'.join(seq_spl)
    name = seq_name + "_" + get_ss_name(c.ss)
    f = open(base_dir + "/" + name + ".pop", "w")
    for i, c in enumerate(clusters):
        mt.add_motif(c.motifs[0])
        node = mt.add_motif(test_helix, end_index=0, end_flip=0)
        if mt.nodes[1].motif.ends[1].r()[1][0] > 0:
            mt.nodes[1].motif.ends[1].flip()
            mt.nodes[1].motif.ends[1].flipped=0
        pop = float(len(c.motifs)) / float(nmotifs)
        energy = -kBT*math.log(pop)
        mt.nodes[1].motif.name = "motif."+str(i)
        f.write("motif."+str(i) + " " + str(len(c.motifs)) + " " + \
                str(pop) + " " + str(energy) + "\n")
        motifs.append(mt.nodes[1].motif)
        mt.remove_node_level()
    f.close()

    mtp = motif_tree_precomputer.MotifTreePrecomputer(name=base_dir+"/"+name,
            max_bps_per_end=0)
    mtp.precompute_motifs(motifs)




def cluster_twoways():
    twoway_lib = motif_library.MotifLibrary(motif_type.TWOWAY)
    twoway_lib.load_all()

    motifs = twoway_lib.motifs()
    start_cluster = SSandSeqCluster(motifs.pop())
    clusters = [ start_cluster ]
    for m in motifs:
        found = 0
        for c in clusters:
            if c.same_seq_and_ss(m):
                c.motifs.append(m)
                found = 1
                break
        if not found:
            clusters.append( SSandSeqCluster(m) )


    mt = motif_tree.MotifTree()
    helix_mlib = motif_library.MotifLibrary(motif_type.HELIX)
    ideal = helix_mlib.get_motif("HELIX.IDEAL")
    origin = target_end_origin(ideal)

    for c in clusters:
        get_two_way_prediction_lib(c, origin, ideal)









if __name__ == '__main__':

    cluster_twoways()
    #get_two_way_prediction_lib()
    exit()
    helix_mlib = motif_library.MotifLibrary(motif_type.HELIX)
    helix_mlib.load_all()
    all_targets = []
    bps = ["AU","UA","CG","GC","GU","UG"]
    seen = []
    #get_bp_step_prediction_lib(["GC","GC"], helix_mlib)
    #exit()
    test_helix = helix_mlib.get_motif("HELIX.IDEAL.2")

    for i,bp1 in enumerate(bps):
        for j,bp2 in enumerate(bps):
            if bp1 + "=" + bp2 in seen:
                continue
            target = [bp1,bp2]
            seen.append(bp1 + "=" + bp2)
            get_bp_step_prediction_lib(target, helix_mlib, test_helix)

