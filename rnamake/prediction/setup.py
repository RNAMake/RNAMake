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


include = """157D
165D
1CSL
1D4R
1DQH
1EHZ
1EVV
1F27
1FUF
1G2J
1I7J
1IK5
1J6S
1J8G
1J9H
1JZV
1KD3
1KD4
1KD5
1KFO
1L2X
1LNT
1MDG
1MSY
1NLC
1NUJ
1NUV
1OSU
1P79
1Q96
1Q9A
1QC0
1QCU
1R3O
1RXB
1T0E
1XPE
1ZCI
1ZEV
255D
259D
2A0P
2A43
2EES
2EET
2EEU
2EEV
2FCX
2FD0
2G32
2G3S
2G91
2G92
2G9C
2GPM
2GQ4
2GQ5
2GQ6
2GQ7
2GRB
2OE5
2OE8
2OEU
2OIY
2PWT
2Q1O
2Q1R
2QEK
2R1S
2R20
2R21
2R22
2V6W
2V7R
2VAL
2VUQ
2W89
2X2Q
2XNZ
2XSL
2ZY6
310D
354D
377D
397D
3BNN
3BNQ
3BNS
3C44
3CGP
3CGS
3CJZ
3CZW
3D0M
3D2V
3DIL
3DS7
3DVV
3DVZ
3DW4
3DW5
3DW6
3DW7
3FO4
3FO6
3GAO
3GER
3GLP
3GM7
3GOT
3GVN
3HGA
3JXQ
3JXR
3LA5
3MEI
3ND3
3ND4
3NJ6
3NJ7
3OK2
3OK4
3P4A
3P4B
3P4C
3P4D
3PDR
3R1D
3R1E
3RG5
3S7C
3S8U
3SD3
3SJ2
3SYW
3TD0
406D
413D
420D
434D
435D
437D
439D
464D
466D
468D
469D
470D
472D
480D
483D
4E58
4E59
4E5C
4E6B
4F8U
4FE5
4FNJ
479D
3MXH
2PN4
4LVZ
3GX5
3RKF
3NKB
1YZD
3CGR
2QUW
3BNO
1QBP
1RNA
2P7D
2Z75
1Y26
4KYY
433D
438D
2HOJ
3B31
3LOA
3SSF
402D
4P95
2GDI
2FQN
2D2L
4K27
405D
3TD1
3E5C
1MWL
421D
3ZP8
3B5S
3TZR
1O9M
1MHK
1GID
398D
3NPQ
4JF2
1Z7F
2AO5
4JAB
3SLQ
280D
3P59
4P5J
3K1V
4ENB
4PQV
3SZX
1DUQ
1T0D
409D
3P22
4K31
353D
4MSR
1I9X
1ZX7""".split("\n")





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


if __name__ == '__main__':
    helix_mlib = motif_library.MotifLibrary(motif_type.HELIX)
    helix_mlib.load_all()
    all_targets = []
    bps = ["AU","UA","CG","GC","GU","UG"]
    seen = []
    #get_bp_step_prediction_lib(["GC","GC"], helix_mlib)
    #exit()
    for i,bp1 in enumerate(bps):
        for j,bp2 in enumerate(bps):
            if bp1 + "=" + bp2 in seen:
                continue
            target = [bp1,bp2]
            seen.append(bp1 + "=" + bp2)
            get_bp_step_prediction_lib(target, helix_mlib)

