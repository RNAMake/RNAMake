import rnamake.motif as motif
import rnamake.motif_tree as motif_tree
import rnamake.util as util
import rnamake.cluster as cluster
import rnamake.settings as settings
import rnamake.setup.motif_library as motif_library
import rnamake.motif_type as motif_type
import rnamake.motif_outputer as motif_outputer
import rnamake.resource_manager as rm
import rnamake.motif_factory as mf
import rnamake.sqlite_library as sqlite_library
import rnamake.motif_ensemble as motif_ensemble
import copy
import math
import itertools

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
    m = mf.factory.motif_from_bps(bps)
    return m

def get_bp_step_prediction_lib_new(helix_mlib):
    not_seen = {}
    mt = motif_tree.MotifTree()
    ideal = helix_mlib.get_motif("HELIX.IDEAL")
    #name = "=".join(targets)
    base_dir = settings.RESOURCES_PATH + "prediction/"

    bps = ["AU","UA","CG","GC","GU","UG"]
    combos = itertools.product(bps, bps)

    mes_keys = ['data', 'name', 'id']
    mes_data = []

    motif_data = []
    motif_keys = ['data', 'name', 'end_name', 'end_id', 'id']

    clusters = []

    for targets in combos:
        aligned_motifs = []
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
                    aligned_motifs.append(m)

        clusters.append([aligned_motifs, targets])

    mes_keys = ['data', 'name', 'id']
    mes_data = []

    motif_data = []
    motif_keys = ['data', 'name', 'end_name', 'end_id', 'id']
    count = 0

    kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
    kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

    for c in clusters:
        motifs = []
        for m in c[0]:
            m_a = mf.factory.can_align_motif_to_end(m, 0)
            if m_a is None:
                continue
            m_a = mf.factory.align_motif_to_common_frame(m_a, 0)
            motifs.append(m_a)

        aligned_clusters = cluster.cluster_motifs(motifs, 0.80)
        target_name = "=".join(c[1])
        clustered_motifs = []
        energies = []

        for j, c_motifs in enumerate(aligned_clusters):
            m = c_motifs.motifs[0]
            m.mtype = motif_type.HELIX
            m.name = target_name + "." + str(j)
            motif_data.append([m.to_str(), m.name, m.ends[0].name(), m.end_ids[0], count])
            count += 1

            pop = float(len(c_motifs.motifs)) / float(len(motifs))
            energy = -kBT*math.log(pop)
            clustered_motifs.append(m)
            energies.append(energy)

        me = motif_ensemble.MotifEnsemble()
        me.setup(clustered_motifs[0].end_ids[0], clustered_motifs, energies)
        mes_data.append([me.to_str(), me.id, count])

        motif = me.members[0].motif
        motif.name =  target_name

        motif_data.append([motif.to_str(), motif.name, motif.ends[0].name(),
                           me.id, count])
        count += 1

    path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/bp_steps_old.db"
    sqlite_library.build_sqlite_library_2(path, mes_data, mes_keys, 'id')
    path = settings.RESOURCES_PATH +"/motif_libraries_new/bp_steps_old.db"
    sqlite_library.build_sqlite_library_2(path, motif_data, motif_keys, 'id')



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
    test_helix = helix_mlib.get_motif("HELIX.IDEAL.12")
    mt.remove_node_level()
    mt.add_motif(test_helix, end_index=1, end_flip=0)
    mt.level += 1
    f = open(base_dir + "/" + name + ".pop", "w")
    lowest = 10000
    lowest_m = None
    for i, c in enumerate(clusters):
        node = mt.add_motif(c.motifs[0], end_index=0)
        if node is None:
            print "made it"
            mt.remove_node_level()
            continue

        node = mt.add_motif(test_helix, end_index=1,  end_flip=0)
        mt.nodes[2].motif.ends[0].flipped=0
        if node is None:
            mt.nodes[2].motif.ends[1].flip()
            mt.nodes[2].motif.ends[1].flipped=0
            node = mt.add_motif(test_helix, end_index=1,  end_flip=0)
        pop = float(len(c.motifs)) / float(nmotifs)
        energy = -kBT*math.log(pop)
        mt.nodes[2].motif.name = name+"."+str(i)
        f.write(name+"."+str(i) + " " + str(len(c.motifs)) + " " + \
                str(pop) + " " + str(energy) + "\n")
        motifs.append(mt.nodes[2].motif)
        motifs[-1]._cache_basepair_frames()
        motifs[-1].structure._cache_coords()
        if lowest > energy:
            lowest = energy
            lowest_m = motifs[-1]
        #motifs[-1].to_pdb('cluster.'+str(i)+'.pdb')
        mt.remove_node_level()
    f.close()

    f = open(base_dir+"/"+name+".new.me", "w")
    i = 0
    for m in motifs:
        m.end_to_add = 0
        m.name += "-0"
        mts = motif_tree_state.motif_to_state_simple(m, 0, 0)
        m2 = motif.str_to_motif(mts.build_string)
        m.to_pdb('cluster.'+str(i)+'.pdb')
        i += 1
        f.write(mts.to_f_str())
    f.close()

    return lowest_m


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

    #mtp = motif_tree_precomputer.MotifTreePrecomputer(name=base_dir+"/"+name,
    #        max_bps_per_end=0)
    #mtp.precompute_motifs(motifs)

    f = open(base_dir+"/"+name+".new.me", "w")
    for m in motifs:
        m.end_to_add = 0
        m.name += "-1"
        mts = motif_tree_state.motif_to_state_simple(m, 0, 0)
        f.write(mts.to_f_str())
    f.close()


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

    #cluster_twoways()
    #get_two_way_prediction_lib()
    #exit()
    helix_mlib = motif_library.MotifLibrary(motif_type.HELIX)
    helix_mlib.load_all()
    all_targets = []
    bps = ["AU","UA","CG","GC","GU","UG"]
    seen = []
    #test_helix = helix_mlib.get_motif("HELIX.IDEAL.2")
    #get_bp_step_prediction_lib(["GC","GC"], helix_mlib)
    get_bp_step_prediction_lib_new(helix_mlib)
    exit()

    ideal_motifs = []

    base_dir = settings.RESOURCES_PATH + "prediction/"

    for i,bp1 in enumerate(bps):
        for j,bp2 in enumerate(bps):
            if bp1 + "=" + bp2 in seen:
                continue
            target = [bp1,bp2]
            name = bp1 + "=" + bp2
            seen.append(name)
            m = get_bp_step_prediction_lib(target, helix_mlib)
            m.name = name
            ideal_motifs.append(m)
    #motif_library_sqlite.build_sqlite_library("bp_steps.db", ideal_motifs)

    f = open(base_dir+"all.new.me", "w")
    for m in ideal_motifs:
        m.end_to_add = 0
        m.name += "-0"
        mts = motif_tree_state.motif_to_state_simple(m, 0, 0)
        m2 = motif.str_to_motif(mts.build_string)
        f.write(mts.to_f_str())
    f.close()



