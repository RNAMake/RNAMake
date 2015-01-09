import rnamake.motif as motif
import rnamake.motif_tree as motif_tree
import rnamake.util as util
import rnamake.cluster as cluster
import rnamake.motif_tree_precomputer as motif_tree_precomputer
import rnamake.settings as settings
import rnamake.motif_library as motif_library
import rnamake.motif_type as motif_type

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
    aligned_motifs = []
    mt = motif_tree.MotifTree()
    ideal = helix_mlib.get_motif("HELIX.IDEAL")
    origin = target_end_origin(ideal)
    name = "=".join(targets)
    base_dir = settings.RESOURCES_PATH + "prediction/"

    for m in helix_mlib.motifs():
        spl = m.name.split(".")
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

    clusters = cluster.cluster_motifs(aligned_motifs)
    motifs = []
    nmotifs = len(aligned_motifs)
    test_helix = helix_mlib.get_motif("HELIX.IDEAL.2")
    f = open(base_dir + "/" + name + ".pop", "w")
    for i, c in enumerate(clusters):
        mt.add_motif(c.motifs[0])
        node = mt.add_motif(test_helix, end_index=0, end_flip=0)
        if node is None:
            mt.nodes[1].motif.ends[1].flip(1)
            mt.nodes[1].motif.ends[1].flipped = 0
        mt.nodes[1].motif.name = name+"."+str(i)
        f.write(name+"."+str(i) + " " + str(len(c.motifs)) + " " + \
                str(float(len(c.motifs))/float(nmotifs))+ "\n")
        motifs.append(mt.nodes[1].motif)
        mt.remove_node_level()
    f.close()
    mtp = motif_tree_precomputer.MotifTreePrecomputer(name=base_dir+"/"+name,
            max_bps_per_end=0)
    mtp.precompute_motifs(motifs)


if __name__ == '__main__':
    helix_mlib = motif_library.MotifLibrary(motif_type.HELIX)
    helix_mlib.load_all()
    all_targets = []
    bps = ["AU","UA","CG","GC","GU","UG"]
    seen = []
    for i,bp1 in enumerate(bps):
        for j,bp2 in enumerate(bps):
            if bp1 + "=" + bp2 in seen:
                continue
            target = [bp1,bp2]
            seen.append(bp1 + "=" + bp2)
            get_bp_step_prediction_lib(target, helix_mlib)

