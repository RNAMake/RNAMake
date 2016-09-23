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
import motif_graph

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
            #if xr is None:
            #    continue
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

    def pose_from_file_new(self, path, gu_are_helix=0, singlet_bp_seperation=0):
        base_motif = motif_factory.factory.motif_from_file(path)
        p = self._copy_info_into_pose(base_motif)
        self._setup_motifs_from_x3dna(p, gu_are_helix, singlet_bp_seperation)
        all_motifs = p.motifs(motif_type.ALL)



        not_helix= []
        helix = []
        for m in all_motifs:
            if m.mtype == motif_type.HELIX:
                helix.append(m)
            else:
                if m.mtype == motif_type.HAIRPIN:
                    m.block_end_add = -1
                not_helix.append(m)

        bp_steps = self._convert_helices_to_bp_steps(helix)
        remove = []
        for bp in bp_steps:
            res = bp.residues()
            for m in not_helix:
                count = 0
                for r in res:
                    if r not in m.residues():
                        break
                    else:
                        count += 1
                if len(res) == count:
                    remove.append(bp)
                    break
        for bp in remove:
            bp_steps.remove(bp)

        all_motifs = not_helix + bp_steps

        #for i, m in enumerate(all_motifs):
        #    m.to_pdb("motif."+str(i)+".pdb")

        start_motif = self._find_best_start_end(p, all_motifs)
        all_motifs.remove(start_motif)
        mg = motif_graph.MotifGraph()
        start_motif.name = p.name + ".motif.0"
        mg.add_motif(start_motif)
        i = 1
        while len(all_motifs) > 0:
            leafs = mg.leafs_and_ends()
            #print len(leafs), len(all_motifs)
            if len(leafs) == 0:
                print "fail"
                mg.write_pdbs()
                exit()

            for l in leafs:
                bp = l[0].data.ends[l[1]]
                found = 0
                next = None
                for m in all_motifs:
                    old_bps = m.get_basepair(name=bp.name())
                    if len(old_bps) == 0:
                        continue
                    pos = m.ends.index(old_bps[0])
                    next = m

                #mg.write_pdbs()
                if next == None:
                    print "fail, next is None"
                    exit()

                next.name = p.name + ".motif." + str(i)
                mg.add_motif(next, parent_index=l[0].index, parent_end_index=l[1])
                all_motifs.remove(next)
                i += 1

        new_p = pose.PoseNew(m=base_motif)
        new_p.mgraph = mg
        return new_p

    def _find_best_start_end(self, p, all_motifs):
        best_m = None
        best_free_ends = []
        for k, m1 in enumerate(all_motifs):
            free_ends = []
            for end in m1.ends:
                free_end_res = []
                not_free_end_res =[]
                for r in end.residues():
                    for c in p.chains():
                        if r == c.first() or r == c.last():
                            free_end_res.append(r)
                    if r not in free_end_res:
                        not_free_end_res.append(r)
                if len(free_end_res) == 2:
                    free_ends.append(end)
                else:
                    #if basepair has single-stranded areas, this will make sure
                    #that this is the first basepair to build from
                    #print len(free_end_res)
                    for r in not_free_end_res:
                        fine = 1
                        cur_chain = None
                        for c in p.chains():
                            if r in c.residues:
                                cur_chain = c
                                break
                        i = cur_chain.residues.index(r)
                        start = 0
                        c_end = i
                        incr = 1
                        if i > len(cur_chain.residues) / 2:
                            start = len(cur_chain.residues)
                            incr = -1
                        while abs(start - c_end) > 1:
                            start += incr
                            new_r = cur_chain.residues[start]
                            bps = p.get_basepair(res1=new_r)
                            for bp in bps:
                                if bp.bp_type == "cW-W":
                                    fine = 0
                                    break
                        if fine:
                            free_end_res.append(r)
                    if len(free_end_res) == 2:
                        free_ends.append(end)
            if best_m is None:
                best_m = m1
                best_free_ends = free_ends
                continue
            if   len(free_ends) > len(best_free_ends):
                best_m = m1
                best_free_ends = free_ends
            elif   (len(best_m.ends) - len(best_free_ends)) > \
                 (len(m1.ends) - len(free_ends)):
                best_m = m1
                best_free_ends = free_ends
            elif (len(best_m.ends) - len(best_free_ends)) > \
                 (len(m1.ends) - len(free_ends)):
                    if len(best_m.ends) > len(m1.ends):
                        best_m = m1
                        best_free_ends = free_ends
        #print best_m, len(best_free_ends)
        return best_m

    def _convert_helices_to_bp_steps(self, helices):
        basepair_steps = []
        for m in helices:
            c = m.chains()[0]
            for i in range(1,len(c.residues)):
                sub_res = [c.residues[i-1], c.residues[i] ]
                sub_bps = []
                for r in sub_res:
                    bp = m.get_basepair(res1=r)[0]
                    for bp_r in bp.residues():
                        if bp_r not in sub_res:
                            sub_res.append(bp_r)
                    if bp not in sub_bps:
                         sub_bps.append(bp)

                bp_step = motif_factory.factory.motif_from_res(sub_res, sub_bps)
                bp_step.mtype = motif_type.HELIX
                basepair_steps.append(bp_step)
        return basepair_steps

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