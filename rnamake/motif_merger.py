import motif_type
import rna_structure
import chain
import secondary_structure_factory as ssf
import secondary_structure
import util

class MotifMerger(rna_structure.RNAStructure):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.residue_map = {}
        self.basepair_map = {}
        self.all_res = {}
        self.all_bp = {}
        self.motifs = {}

    def copy(self, new_motifs=None):
        new_all_res = {}
        new_all_bp = {}
        motifs = {}
        for m in new_motifs:
            motifs[m.id] = m
            for r in m.residues():
                new_all_res[r.uuid] = r
            for bp in m.basepairs:
                new_all_bp[bp.uuid] = bp

        new_struct = MotifMerger()
        new_struct.residue_map = self.residue_map

        chains = []
        for c in self.chains():
            new_res = []
            for r in c.residues:
                new_res.append(new_all_res[r.uuid])
            chains.append(chain.Chain(new_res))

        new_bps = []
        for bp in self.basepairs:
            new_bp = new_all_bp[bp.uuid]
            new_bps.append(new_bp)

        new_ends = []
        for end in self.ends:
            new_end = new_all_bp[end.uuid]
            new_ends.append(new_end)

        new_struct.structure.chains = chains
        new_struct.basepairs = new_bps
        new_struct.ends = new_ends
        new_struct.all_res = new_all_res
        new_struct.all_bp = new_all_bp
        new_struct.motifs = motifs

        return new_struct

    def add_motif(self, m, m_end=None, parent=None, parent_end=None):
        for r in m.residues():
            self.all_res[r.uuid] = r
        for bp in m.basepairs:
            self.all_bp[bp.uuid] = bp
        self.motifs[m.id] = m

        if parent is None:
           self.structure.chains.extend([c.subchain(0) for c in m.chains() ])
           self.basepairs.extend(m.basepairs)
           self.ends.extend(m.ends)

        else:
            self._merge_motifs(m, parent, m_end, parent_end,
                               [c.subchain(0) for c in m.chains()],
                               self.chains())


    def connect_motifs(self, m1, m2, m1_end, m2_end):
        self._merge_motifs(m1, m2, m1_end, m2_end,self.chains(), self.chains())

    def remove_motif(self, m):
        removed_res = self._remove_res_from_chains(m)
        keep_bps = []
        removed_bp = []
        for bp in self.basepairs:
            if bp.res1 in removed_res and bp.res2 in removed_res:
                removed_bp.append(bp)
                continue
            keep_bps.append(bp)
        self.basepairs = keep_bps
        self.ends = rna_structure.ends_from_basepairs(self.structure, self.basepairs)
        self.motifs.pop(m.id, None)

        #remove stuff from dictionaries
        for r in removed_res:
            self.all_res.pop(r.uuid, None)
            if r.uuid in self.residue_map:
                self.residue_map.pop(r.uuid, None)
            key = self._get_key_for_value(r.uuid)
            if key is not None:
                self.residue_map.pop(key, None)

        for bp in removed_bp:
            self.all_bp.pop(bp, None)

    def update_motif(self, m):
        for r in m.residues():
            act_r = self.get_residue(uuid=r.uuid)
            if r.uuid != act_r.uuid:
                continue

            for c in self.structure.chains:
                index = -1
                for i, r1 in enumerate(c.residues):
                    if r.uuid == r1.uuid and r.name != 'N':
                        index = i
                        break

                if index != -1:
                    c.residues[index] = r
        self.motifs[m.id] = m

    def secondary_structure(self):

        ss = ssf.factory.secondary_structure_from_motif(self)
        ss_motifs = []

        for m in self.motifs.values():
            ss_chains = []
            for c in m.chains():
                ss_res = []
                for r in c.residues:
                    correct_r = self.get_residue(uuid=r.uuid)
                    ss_r = ss.get_residue(uuid=correct_r.uuid)
                    ss_res.append(ss_r)
                ss_chains.append(secondary_structure.Chain(ss_res))
            ss_struct = secondary_structure.Structure(ss_chains)
            ss_bps = []
            for bp in m.basepairs:
                if bp.bp_type != "cW-W":
                    continue
                if not util.wc_bp(bp) and not util.gu_bp(bp):
                    continue
                correct_bp = self.get_basepair(bp_uuid=bp.uuid)
                ss_bp = ss.get_bp(uuid=correct_bp[0].uuid)

                if ss_bp is None:
                    print correct_bp[0].uuid
                    print correct_bp[0], ss_bp, "type: ", correct_bp[0].bp_type
                    print correct_bp[0].res1 in self.residues(), correct_bp[0].res2 in self.residues()
                    #raise ValueError("did not find a bp we should of")

            ss_rna_struct = secondary_structure.RNAStructure(
                                ss_struct, ss_bps, [], m.name, m.path,
                                m.mtype, m.score, m.end_ids)

            ss_ends = []
            for end in m.ends:
                correct_bp = self.get_basepair(bp_uuid=end.uuid)
                ss_bp = ss.get_bp(uuid=correct_bp[0].uuid)
                ss_ends.append(ss_bp)
            ss_rna_struct.ends = ss_ends
            ss_motifs.append(secondary_structure.Motif(r_struct=ss_rna_struct,
                                                       id=m.id))

        ss_struct = secondary_structure.Structure(ss.chains)
        ss_p = secondary_structure.Pose(ss_struct, ss.basepairs, ss.ends)
        ss_p.residue_map = self.residue_map
        ss_p.motifs = ss_motifs

        return ss_p

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        """
        wrapper to self.structure.get_residue()
        """
        r = super(self.__class__, self).get_residue(num, chain_id, i_code, uuid)

        if r is not None:
            return r
        if uuid is not None and uuid in self.residue_map:
            return super(self.__class__, self).get_residue(num, chain_id, i_code,
                                                           self.residue_map[uuid])

    def get_basepair(self, bp_uuid=None, res1=None, res2=None, uuid1=None,
                     uuid2=None, name=None):
        bps = super(self.__class__, self).get_basepair(bp_uuid, res1, res2, uuid1, uuid2, name)
        if len(bps) == 0:
            if bp_uuid in self.basepair_map:
                new_uuid = self.basepair_map[bp_uuid]
                return super(self.__class__, self).get_basepair(bp_uuid=new_uuid)
            else:
                return bps

        else:
            return bps

    def _remove_res_from_chains(self, m):
        chains = self.chains()
        new_chains = []
        removed_res = []
        removed = 1
        count = 0
        while removed:
            removed = 0
            split = []
            for c in chains:
                for c1 in m.chains():
                    found = []
                    start = -1
                    for r in c1.residues:
                        if r in c.residues:
                            if start == -1:
                                start = c.residues.index(r)
                            found.append(r)

                    if start == -1:
                        continue

                    next_chains = []
                    if start != 0:
                        new_c = c.subchain(0, start)
                        if found[0].uuid in self.residue_map.values():
                            other_uuid = self._get_key_for_value(found[0].uuid)
                            new_r = self.all_res[other_uuid]
                            new_c.residues.insert(0, new_r)
                        next_chains.append(new_c)

                    if start+len(found) != len(c):
                        new_c = c.subchain(start+len(found))
                        new_chains.append(new_c)
                        if found[-1].uuid in self.residue_map.values():
                            other_uuid = self._get_key_for_value(found[-1].uuid)
                            new_r = self.all_res[other_uuid]
                            new_c.residues.append(new_r)

                    if start != 0 or len(c) != start+len(found):
                        split.append(c)

                    if len(next_chains) > 0:
                        new_chains.extend(next_chains)
                    removed_res.extend(found)
                    removed = 1

            if len(new_chains) == 0 and count == 0:
                chains = []
                break

            if len(new_chains) == 0:
                break

            #unused chains to still be included
            for c in chains:
                if c not in split:
                    new_chains.append(c)

            chains = new_chains
            new_chains = []
            count += 1

        self.structure.chains = chains
        res = self.residues()
        for bp in self.all_bp.values():
            if bp in self.basepairs:
                continue
            if bp.res1 in res and bp.res2 in res:
                self.basepairs.append(bp)

        return removed_res

    def _get_key_for_value(self, uuid):
        for k,v in self.residue_map.iteritems():
            if v == uuid:
                return k
        return None

    def _merge_motifs(self, m1, m2, m1_end, m2_end, m1_chains, m2_chains):
        s_chain_map = rna_structure.get_chain_end_map(m2_chains, m2_end)
        m_chain_map = rna_structure.get_chain_end_map(m1_chains, m1_end)

        if   m2.mtype == motif_type.HELIX and m1.mtype != motif_type.HELIX:
            merged_chains = self._merge_chains(m_chain_map, s_chain_map)
            self.basepair_map[m2_end.uuid] = m1_end.uuid
            self.basepairs.remove(m2_end)
        else:
            merged_chains = self._merge_chains(s_chain_map, m_chain_map)
            self.basepair_map[m1_end.uuid] = m2_end.uuid

        for c in s_chain_map.chains() + m_chain_map.chains():
            if c in self.structure.chains:
                self.structure.chains.remove(c)
        for c in merged_chains:
            if c is not None:
                self.structure.chains.append(c)
        for c in m1_chains:
            if c not in self.structure.chains and c not in m_chain_map.chains():
                self.structure.chains.append(c)

        for e in m1.ends:
            if e != m1_end and e not in self.ends:
                self.ends.append(e)

        #update ends
        self.ends.remove(m2_end)
        try:
            self.ends.remove(m1_end)
        except:
            pass

        #update basepairs
        uuids = [r.uuid for r in self.residues()]
        for bp in m1.basepairs:
            if bp.res1.uuid in uuids and bp.res2.uuid in uuids and \
               bp not in self.basepairs:
                self.basepairs.append(bp)

        if m1_chains == m2_chains:
            self.basepairs.remove(m1_end)


    def _merge_chains(self, cm1, cm2):
        merged_chain_1, merged_chain_2 = None, None
        if   cm1.is_hairpin() and cm2.is_hairpin():
            raise ValueError("cannot merge an hairpin with another hairpin")
        elif cm1.is_hairpin():
            merged_chain_1 = self._get_merged_hairpin(cm2.p3_chain, cm2.p5_chain,
                                                      cm1.p5_chain, 0, 1)
        elif cm2.is_hairpin():
            merged_chain_1 = self._get_merged_hairpin(cm1.p5_chain, cm1.p3_chain,
                                                      cm2.p5_chain, 1, 1)
        else:
            merged_chain_1 = self._get_merged_chain(cm1.p5_chain, cm2.p3_chain, 1, 1)
            merged_chain_2 = self._get_merged_chain(cm1.p3_chain, cm2.p5_chain, 0, 1)

        return merged_chain_1, merged_chain_2

    def _get_merged_hairpin(self, c1, c2, hairpin, join_by_3prime=0,
                            remove_overlap=0):
        merged_chain = self._get_merged_chain(c1, hairpin, join_by_3prime,
                                              remove_overlap)
        merged_chain = self._get_merged_chain(merged_chain, c2, join_by_3prime,
                                              remove_overlap)
        return merged_chain

    def _get_merged_chain(self, c1, c2, join_by_3prime=0, remove_overlap=0):
        """
        Merges two chains together that share a common resiude

        :param c1: chain 1
        :param c2: chain 2
        :param join_by_3_prime: joins in the 3prime direction instead of the
            standard (optional)
        :param remove_overlap: removes the overlap residue between the two
            chains (optional)
        :type c1: chain object

        Example:

        chain1            chain2
        5'_|_|_|_|_|_3' + 5'_|_|_|_|_|_3' =
        5'_|_|_|_|_|_|_|_|_|_|_ 3'

        Notice one residue is lost at 5' end of the chain2, this is the overlap
        residue which was used to align the chains during alignment, you can
        set remove_overlap to 0 to stop that from happning

        """
        merged_chain = chain.Chain()
        chain1_res, chain2_res = c1.residues, c2.residues
        if join_by_3prime:
            chain1_res, chain2_res = chain1_res[::-1], chain2_res[::-1]
        merged_chain.residues = list(chain1_res)
        if remove_overlap:
            r = chain2_res.pop(0)
            self.residue_map[r.uuid] = chain1_res[-1].uuid
        merged_chain.residues.extend(list(chain2_res))
        if join_by_3prime:
            merged_chain.residues = merged_chain.residues[::-1]
        return merged_chain
