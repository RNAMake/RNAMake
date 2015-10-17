import motif_factory
import secondary_structure_factory
import secondary_structure
import motif_type
import pose
import x3dna
import util
import chain
import motif
import structure

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

    def pose_from_motif_tree_old(self, structure, basepairs, motifs, designable):

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


        #update all motifs in pose, make sure that the correct basepairs are used
        #since basepairs are used to align there are two choices of which basepair to
        #be used
        for j, m in enumerate(motifs):
            bps = []
            residue_map = {}
            for bp in m.basepairs:
                new_bp = p.get_basepair(bp_uuid=bp.uuid)
                #print bp.residues(), len(new_bp)
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
            #print m.name
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
                        raise ValueError("cannot find residue")
                chains.append(chain.Chain(res))

            if len(bps) != len(m.basepairs):
                raise ValueError("something went horribly wrong: did not find all basepairs")

            m_copy = motif.Motif()
            m_copy.mtype = m.mtype
            m_copy.name = m.name
            m_copy.structure.chains = chains
            m_copy.basepairs = bps
            m_copy.ends = motif_factory.factory._setup_basepair_ends(m_copy.structure, bps)
            motif_factory.factory._setup_secondary_structure(m_copy)
            if m_copy.mtype is not motif_type.HELIX:
                best_ends = []
                best_end_ids = []
                for i, end in enumerate(m.ends):
                    best = 1000
                    best_end = None
                    best_end_id = None
                    for c_end in m_copy.ends:
                        dist = util.distance(end.d(), c_end.d())
                        if dist < best:
                            best_end = c_end
                            best_end_id = m.end_ids[i]
                            best = dist
                    best_ends.append(best_end)
                    best_end_ids.append(best_end_id)

                m_copy.ends = best_ends
                m_copy.end_ids = best_end_ids


            p.motif_dict[m_copy.mtype].append(m_copy)
            p.motif_dict[motif_type.ALL].append(m_copy)

        self._add_secondary_structure_motifs(p)
        self.standardize_pose(p)
        return p

    def pose_from_motif_tree(self, structure, basepairs, nodes, designable):

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


        #update all motifs in pose, make sure that the correct basepairs are used
        #since basepairs are used to align there are two choices of which basepair to
        #be used
        for j, n in enumerate(nodes):
            m = n.data
            bps = []
            residue_map = {}
            for bp in m.basepairs:
                new_bp = p.get_basepair(bp_uuid=bp.uuid)
                if len(new_bp) != 0:
                    bps.append(new_bp[0])
                    continue
                best = 1000
                best_bp = None
                for c in n.connections:
                    if c is None:
                        continue
                    if bp not in m.ends:
                        raise ValueError("cannot find a non end basepair")
                    if c.end_index(n.index) != m.ends.index(bp):
                        continue
                    partner = c.partner(n.index)
                    alt_bp = partner.data.ends[c.end_index(partner.index)]
                    pose_bps = p.get_basepair(bp_uuid=alt_bp.uuid)
                    if len(pose_bps) == 0:
                        raise ValueError("cannot find partner end bp in pose")
                    #TODO Check if this is okay, not sure if I have to try and match
                    #TODO residues correctly
                    residue_map[bp.res1] = pose_bps[0].res1
                    residue_map[bp.res2] = pose_bps[0].res2
                    bps.append(pose_bps[0])
                    break

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
                        raise ValueError("cannot find residue")
                chains.append(chain.Chain(res))
            if len(bps) != len(m.basepairs):
                raise ValueError("something went horribly wrong: did not find all basepairs")

            m_copy = motif.Motif()
            m_copy.mtype = m.mtype
            m_copy.name = m.name
            m_copy.structure.chains = chains
            m_copy.basepairs = bps
            m_copy.ends = motif_factory.factory._setup_basepair_ends(m_copy.structure, bps)
            motif_factory.factory._setup_secondary_structure(m_copy)
            if m_copy.mtype is not motif_type.HELIX:
                best_ends = []
                best_end_ids = []
                for i, end in enumerate(m.ends):
                    best = 1000
                    best_end = None
                    best_end_id = None
                    for c_end in m_copy.ends:
                        dist = util.distance(end.d(), c_end.d())
                        if dist < best:
                            best_end = c_end
                            best_end_id = m.end_ids[i]
                            best = dist
                    best_ends.append(best_end)
                    best_end_ids.append(best_end_id)

                m_copy.ends = best_ends
                m_copy.end_ids = best_end_ids


            p.motif_dict[m_copy.mtype].append(m_copy)
            p.motif_dict[motif_type.ALL].append(m_copy)

        self._add_secondary_structure_motifs(p)
        self.standardize_pose(p)
        return p

    def pose_from_motif_tree_new(self, structure, basepairs, nodes, designable):
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

        #update all motifs in pose, make sure that the correct basepairs are used
        #since basepairs are used to align there are two choices of which basepair to
        #be used
        p_res = p.residues()
        for j, n in enumerate(nodes):
            m = n.data
            residue_map = {}
            print j
            chains = []
            res = []
            for c in m.chains():
                bounds = [0, len(c.residues)]
                first = 0
                start_pos = 0
                found = 0
                for i, r in enumerate(c.residues):
                    r_new = p.get_residue(uuid=r.uuid)
                    if r_new is not None:
                        first = i
                        start_pos = p_res.index(r_new)
                        found = 1
                        break

                bounds[0] = start_pos - first
                bounds[1] = bounds[0] + len(c.residues)
                print bounds[0], bounds[1]

                chains.append(chain.Chain(p_res[bounds[0]:bounds[1]]))
                res.extend(chains[-1].residues)

            bps = self._subselect_bps_from_res(p.basepairs, res)
            m_copy = motif_factory.factory.motif_from_chains(chains, bps)
            m_copy.to_pdb("m."+str(j)+".pdb")
            m.to_pdb("org."+str(j)+".pdb")
            continue

            exit()

            if len(bps) != len(m.basepairs):
                raise ValueError("something went horribly wrong: did not find all basepairs")

            m_copy = motif.Motif()
            m_copy.mtype = m.mtype
            m_copy.name = m.name
            if m_copy.mtype is not motif_type.HELIX:
                best_ends = []
                best_end_ids = []
                for i, end in enumerate(m.ends):
                    best = 1000
                    best_end = None
                    best_end_id = None
                    for c_end in m_copy.ends:
                        dist = util.distance(end.d(), c_end.d())
                        if dist < best:
                            best_end = c_end
                            best_end_id = m.end_ids[i]
                            best = dist
                    best_ends.append(best_end)
                    best_end_ids.append(best_end_id)

                m_copy.ends = best_ends
                m_copy.end_ids = best_end_ids


            p.motif_dict[m_copy.mtype].append(m_copy)
            p.motif_dict[motif_type.ALL].append(m_copy)

        self._add_secondary_structure_motifs(p)
        self.standardize_pose(p)

        return p

    def pose_from_motif_graph(self):
        pass


    def _old_pose_from_motif_tree(self):

        return

        """
        for bp in m.basepairs:
            new_bp = p.get_basepair(bp_uuid=bp.uuid)
            if len(new_bp) != 0:
                bps.append(new_bp[0])
                continue
            best = 1000
            best_bp = None
            for c in n.connections:
                if c is None:
                    continue
                if bp not in m.ends:
                    raise ValueError("cannot find a non end basepair")
                if c.end_index(n.index) != m.ends.index(bp):
                    continue
                partner = c.partner(n.index)
                alt_bp = partner.data.ends[c.end_index(partner.index)]
                pose_bps = p.get_basepair(bp_uuid=alt_bp.uuid)
                if len(pose_bps) == 0:
                    raise ValueError("cannot find partner end bp in pose")
                #TODO Check if this is okay, not sure if I have to try and match
                #TODO residues correctly
                residue_map[bp.res1] = pose_bps[0].res1
                residue_map[bp.res2] = pose_bps[0].res2
                bps.append(pose_bps[0])
                break

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
                    raise ValueError("cannot find residue")
            chains.append(chain.Chain(res))
        """

    def _subselect_bps_from_res(self, bps, res):
        sub_bps = []

        uuid_dict = { r.uuid : 1 for r in res }

        for bp in bps:
            if bp.res1.uuid in uuid_dict and bp.res2.uuid in uuid_dict:
                sub_bps.append(bp)

        return sub_bps

    def split_pose_by_chains(self, p):
        ps = []
        for c in p.chains():
            p_new = pose.Pose()
            s = structure.Structure([c.copy()])
            p.structure = s
            bps = []
            for bp in p.basepairs:
                pass

            ps.append(p)

        return ps

    #TODO try and minimize clashes
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