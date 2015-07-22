import motif_factory
import secondary_structure_factory
import secondary_structure
import motif_type
import pose
import x3dna
import util
import chain
import motif

class PoseFactory(object):
    def __init__(self):
        pass

    def _copy_info_into_pose(self, m):
        p = pose.Pose()
        p.structure = m.structure
        p.ends      = m.ends
        p.end_ids   = m.end_ids
        p.basepairs = m.basepairs
        p.name      = m.name
        p.path      = m.path
        p.beads     = m.beads
        p.score     = m.score
        p.secondary_structure = m.secondary_structure

        return p

    #TODO implement finding tertiary contacts
    def _setup_motifs_from_x3dna(self, p, gu_are_helix, singlet_bp_seperation):
        x = x3dna.X3dna()
        x3dna_motifs = x.get_motifs(p.path)
        motifs = []
        for xm in x3dna_motifs:
            m = self._convert_x3dna_to_motif(xm, p)
            m.mtype = xm.mtype
            motifs.append(m)

        do_not_keep = {}
        pairs = []
        for i, m in enumerate(motifs):
            keep = 0

            if gu_are_helix and m.mtype == motif_type.TWOWAY:
                keep = self._remove_gu_basepair_motifs(m)
                if not keep:
                    do_not_keep[m] = 1

            if not singlet_bp_seperation and m.mtype == motif_type.TWOWAY:
                pair = self._merge_motifs_with_single_bp(m, motifs, i)
                if pair is not None:
                    do_not_keep[pair[0]] = 1
                    do_not_keep[pair[1]] = 1
                    pairs.append(pair)

        for m in motifs:
            if m not in do_not_keep:
                p.motif_dict[m.mtype].append(m)
                p.motif_dict[motif_type.ALL].append(m)

        for pair in pairs:
            res = []
            basepairs = []
            res.extend(pair[0].residues())
            basepairs.extend(pair[0].basepairs)
            for r in pair[1].residues():
                if r not in res:
                    res.append(r)
            for bp in pair[1].basepairs:
                if bp not in basepairs:
                    basepairs.append(bp)
            m = motif_factory.factory.motif_from_res(res, basepairs)
            m.mtype = motif_type.TWOWAY
            p.motif_dict[m.mtype].append(m)
            p.motif_dict[motif_type.ALL].append(m)

    def _remove_gu_basepair_motifs(self, m):
        res = []
        keep = 0
        for bp in m.basepairs:
            for r in bp.residues():
                if r not in res:
                    res.append(r)
            if not (util.wc_bp(bp) and bp.bp_type == "cW-W" or util.gu_bp(bp)):
                keep = 1

        if len(res) != len(m.residues()):
            keep = 1
        return keep

    def _merge_motifs_with_single_bp(self, m, motifs, i):
        for j, m2 in enumerate(motifs):
            if i >= j or m2.mtype != motif_type.TWOWAY:
                continue
            paired = 0
            for end in m.ends:
                if end in m2.ends:
                    return [m, m2]

            for c1 in m.chains():
                for c2 in m2.chains():
                    if abs(c1.first().num -c2.last().num) == 1:
                        return [m, m2]
                    if abs(c2.first().num -c1.last().num) == 1:
                        return [m, m2]

        return None

    def _convert_x3dna_to_motif(self, xm, p):
        res = []
        for xr in xm.residues:
            r = p.get_residue(num=xr.num, chain_id=xr.chain_id, i_code=xr.i_code)
            res.append(r)
        basepairs = []
        for r in res:
            bps = p.get_basepair(res1=r)
            for bp in bps:
                if bp.res1 in res and bp.res2 in res and bp not in basepairs:
                    basepairs.append(bp)

        m = motif_factory.factory.motif_from_res(res, basepairs)
        return m

    def pose_from_file(self, path, gu_are_helix=1, singlet_bp_seperation=0):
        base_motif = motif_factory.factory.motif_from_file(path)
        p = self._copy_info_into_pose(base_motif)
        self._setup_motifs_from_x3dna(p, gu_are_helix, singlet_bp_seperation)

        return p

    def _add_secondary_structure_motifs(self, p):
        ss = p.secondary_structure
        for m in p.motifs(motif_type.ALL):
            ss_ends, ss_bps, ss_chains = [], [], []
            for c in m.chains():
                ss_res = []
                for r in c.residues:
                    ss_r = ss.get_residue(uuid=r.uuid)
                    ss_res.append(ss_r)
                ss_chains.append(secondary_structure.Chain(ss_res))
            for bp in m.basepairs:
                res1 = ss.get_residue(uuid=bp.res1.uuid)
                res2 = ss.get_residue(uuid=bp.res2.uuid)
                ss_bp = ss.get_bp(res1, res2)
                if ss_bp is None:
                    continue
                ss_bps.append(ss_bp)
            ss_end_ids = []
            for i, end in enumerate(m.ends):
                res1 = ss.get_residue(uuid=end.res1.uuid)
                res2 = ss.get_residue(uuid=end.res2.uuid)
                new_bp = ss.get_bp(res1, res2)
                if new_bp is None:
                    raise ValueError("cannot find ss end")
                ss_ends.append(new_bp)
                ss_end_ids.append(m.end_ids[i])

            type = motif_type.type_to_str(m.mtype)
            ss_motif = secondary_structure.SecondaryStructureMotif(type, ss_ends, ss_chains)
            ss_motif.basepairs = ss_bps
            if type not in ss.elements:
                ss.elements[type] = []
            ss_motif.name = m.name
            ss_motif.end_ids = ss_end_ids
            ss.elements[type].append(ss_motif)
            ss.elements['ALL'].append(ss_motif)

    def pose_from_motif_tree(self, structure, basepairs, motifs, designable):

        p = pose.Pose()
        p.designable = designable
        p.name = "assembled"
        p.path = "assembled"
        p.structure = structure
        p.basepairs = basepairs
        p.ends = motif_factory.factory._setup_basepair_ends(structure, basepairs)
        ss = secondary_structure_factory.factory.secondary_structure_from_motif(p)
        p.secondary_structure = ss
        p.end_ids = ["" for x in p.ends]
        for i, end in enumerate(p.ends):
            res1 = ss.get_residue(uuid=end.res1.uuid)
            res2 = ss.get_residue(uuid=end.res2.uuid)
            ss_end = ss.get_bp(res1, res2)
            p.end_ids[i] = secondary_structure.assign_end_id(ss, ss_end)

        for m in motifs:
            bps = []
            residue_map = {}
            for bp in m.basepairs:
                new_bp = p.get_basepair(bp_uuid=bp.uuid)
                if len(new_bp) != 0:
                    bps.append(new_bp[0])
                    continue

                best = 1000
                best_bp = None
                for m2 in motifs:
                    if m == m2:
                        continue
                    for end in m2.ends:
                        alt_bp = p.get_basepair(bp_uuid=end.uuid)
                        if len(alt_bp) == 0:
                            continue
                        dist = util.distance(bp.d(), alt_bp[0].d())
                        if dist < best:
                            best_bp = alt_bp[0]
                            best = dist
                if best_bp is None:
                    continue
                new_bp = [best_bp]
                for r1 in best_bp.residues():
                    r1_cent = util.center(r1.atoms)
                    best_match = None
                    best_dist = 1000
                    for r2 in bp.residues():
                        r2_cent = util.center(r2.atoms)
                        dist = util.distance(r1_cent, r2_cent)
                        if dist < best_dist:
                            best_dist = dist
                            best_match = r2
                    residue_map[best_match] = r1
                bps.append(new_bp[0])

            chains = []
            for c in m.chains():
                res = []
                for r in c.residues:
                    new_r = p.get_residue(uuid=r.uuid)
                    if new_r is not None:
                        res.append(new_r)
                    elif r in residue_map:
                        res.append(residue_map[r])
                    else:
                        print r, r.uuid
                        #raise ValueError("cannot find residue")
                chains.append(chain.Chain(res))

            m_copy = motif.Motif()
            m_copy.mtype = m.mtype
            m_copy.name = m.name
            m_copy.structure.chains = chains
            m_copy.basepairs = bps
            m_copy.ends = motif_factory.factory._setup_basepair_ends(m_copy.structure, bps)
            motif_factory.factory._setup_secondary_structure(m_copy)
            if m_copy.mtype is not motif_type.HELIX:
                m_copy.end_ids = m.end_ids
            p.motif_dict[m_copy.mtype].append(m_copy)
            p.motif_dict[motif_type.ALL].append(m_copy)

        self._add_secondary_structure_motifs(p)
        self.standardize_pose(p)
        return p

    def standardize_pose(self, p):
        fail = 0
        added_helix = motif_factory.factory.added_helix
        for i in range(len(p.ends)):
            m_added = motif.get_aligned_motif(p.ends[i],
                                              added_helix.ends[0],
                                              added_helix)

            if not motif_factory.factory._steric_clash(p, m_added):
                continue

            p.ends[i].flip()

            m_added = motif.get_aligned_motif(p.ends[i],
                                              added_helix.ends[0],
                                              added_helix)



factory = PoseFactory()