import x3dna
import motif_type
import structure
import chain
import secondary_structure
import motif_state
import util
import graph
import user_warnings
from primitives.rna_structure import ends_from_basepairs, assign_end_id, end_id_to_seq_and_db


class ChainNodeData(object):
    def __init__(self, c, m_id, prime5_override=0, prime3_override=0):
        self.c = c
        self.m_id = m_id
        self.prime5_override = prime5_override
        self.prime3_override = prime3_override

    def included_res(self):
        res = [ r for r in  self.c]
        if self.prime5_override:
            res.pop(0)
        if self.prime3_override:
            res.pop()
        return res


class MotifMergerType(object):
    # sequence matters for these motifs
    SPECIFIC_SEQUENCE = 0
    # sequence does not matter such as idealized helices
    NON_SPECIFIC_SEQUENCE = 1


class Merger(object):
    __slots__ = [
        "_all_bps",
        "_motifs",
        "_res_overrides",
        "_bp_overrides",
        "_rebuild_structure",
        "_rna_structure",
        "_chain_graph",
        "_aligned",
        "_m_pos",
    ]

    def __init__(self):
        self._all_bps = {}
        self._motifs = {}
        self._res_overrides = {}
        self._bp_overrides = {}

        self._aligned = {}
        self._m_pos = {}
        self._rebuild_structure = 1
        self._chain_graph = graph.GraphStatic()
        self._rna_structure = None

    def add_motif(self, m, m_end=None, parent=None, parent_end=None):
        new_chains = [c for c in m.get_chains()]

        for c in new_chains:
            data = ChainNodeData(c, m.uuid)
            self._chain_graph.add_data(data, n_children=2, orphan=1)

        for bp in m.iter_basepairs():
            self._all_bps[bp.uuid] = bp
        for end in m.iter_ends():
            self._all_bps[end.uuid] = end

        if parent is not None:
            self._aligned[m.uuid] = 1
            self._link_motifs(parent, m, parent_end, m_end)
        else:
            self._aligned[m.uuid] = 0

        self._m_pos[m.uuid] = len(self._motifs)
        self._motifs[m.uuid] = m
        self._rebuild_structure = 1

    def remove_motif(self, m):
        for end in m.ends:
            if end.uuid in self.bp_overrides:
                del self.bp_overrides[end.uuid]
            remove = []
            for uuid1, uuid2 in self.bp_overrides.iteritems():
                if uuid2 == end.uuid:
                    remove.append(uuid1)
            for r in remove:
                self.bp_overrides.pop(r, None)

        for bp in m.basepairs:
            if bp.uuid in self.all_bps:
                self.all_bps.pop(bp.uuid, None)

        remove = []
        for n in self.chain_graph.nodes:
            if n.data.m_id == m.id:
                remove.append(n)

        for r in remove:
            for c in r.connections:
                if c is None:
                    continue
                p = c.partner(r.index)
                p_i = c.end_index(p.index)

                if   p_i == 0 and p.data.prime5_override == 1:
                    del self.res_overrides[p.data.c.first().uuid]
                    p.data.prime5_override = 0
                elif p_i == 1 and p.data.prime3_override == 1:
                    del self.res_overrides[p.data.c.last().uuid]
                    p.data.prime3_override = 0
            self.chain_graph.remove_node(r.index)

        del self.motifs[m.id]
        self.rebuild_structure = 1

    def connect_motifs(self, m1, m2, m1_end, m2_end):
        self._link_motifs(m1, m2, m1_end, m2_end)

    def __get_ss_structure(self, m, ss):
        ss_chains = []
        for c in m.get_chains():
            ss_res = []
            for r in c:
                r_cur = r
                if r.uuid in self._res_overrides:
                    r_cur = self.get_residue(self._res_overrides[r.uuid])
                ss_r = ss.get_residue(uuid=r_cur.uuid)
                if ss_r is None:
                    raise ValueError("could not find residue during ss build: " +
                                        str(r_cur) + "from " + m.name)
                ss_res.append(ss_r)
            ss_chains.append(secondary_structure.Chain(ss_res))

        res = []
        chain_cuts = []
        for i, c in enumerate(ss_chains):
            for r in c:
                res.append(r)
            if i < len(ss_chains) - 1:
                chain_cuts.append(len(res))
        chain_cuts.append(len(res))

        ss_struct = secondary_structure.Structure(res, chain_cuts)
        return ss_struct

    def __get_ss_basepairs(self, m, ss):
        ss_bps = []
        for bp in m.iter_basepairs():
            if bp.bp_type != x3dna.X3dnaBPType.cWUW:
                continue
            if not util.wc_bp(bp, m) and not util.gu_bp(bp, m):
                continue
            correct_bp = bp
            if bp.uuid in self._bp_overrides:
                correct_bp = self.get_basepair(self._bp_overrides[bp.uuid])
            ss_bp = ss.get_basepair(bp_uuid=correct_bp.uuid)
            if ss_bp is None:
                raise ValueError("could not find basepair during ss build")
            ss_bps.append(ss_bp)
        return ss_bps

    def __get_ss_ends(self, m, ss):
        ss_ends = []
        for end in m.iter_ends():
            correct_bp = end
            if end.uuid in self._bp_overrides:
                correct_bp = self.get_basepair(self._bp_overrides[end.uuid])
            ss_bp = ss.get_end(bp_uuid=correct_bp.uuid)
            if ss_bp is None:
                ss_bp = ss.get_basepair(bp_uuid=correct_bp.uuid)
            if ss_bp is None:
                raise ValueError("cannot not find end durign ss build")
            ss_ends.append(ss_bp)
        return ss_ends

    def get_merged_secondary_structure(self):
        ss = self.get_merged_structure().get_secondary_structure()

        ss_motifs = []
        for m in self._motifs.values():
            ss_struct = self.__get_ss_structure(m, ss)
            ss_bps = self.__get_ss_basepairs(m, ss)
            ss_ends = self.__get_ss_ends(m, ss)
            end_ids = [ assign_end_id(ss_struct, ss_bps, ss_ends, end) for end in ss_ends]

            ss_motif = secondary_structure.Motif(ss_struct, ss_bps, ss_ends,
                                                 end_ids, m.mtype, m.uuid)
            ss_motifs.append(ss_motif)

        ss_chains = ss.get_chains()

        res = []
        chain_cuts = []
        for i, c in enumerate(ss_chains):
            for r in c:
                res.append(r)
            if i < len(ss_chains) - 1:
                chain_cuts.append(len(res))
        chain_cuts.append(len(res))

        ss_struct  = secondary_structure.Structure(res, chain_cuts)
        ss_bps     = [ bp for bp in ss.iter_basepairs()]
        ss_ends    = [ end for end in ss.iter_ends()]
        ss_end_ids = [ ss.get_end_id(i) for i in range(ss.num_ends())]

        return secondary_structure.Pose(ss_struct, ss_bps, ss_ends, ss_end_ids,
                                        ss_motifs)

    def get_merged_structure(self):
        if self._rebuild_structure == 1:
            self._build_structure()
            self._rebuild_structure = 0

        return self._rna_structure

    def get_residue(self, uuid):
        for m in self._motifs.values():
            r = m.get_residue(uuid=uuid)
            if r is not None:
                return r
        return None

    def get_basepair(self, uuid):
        if uuid in self._all_bps:
            return self._all_bps[uuid]
        else:
            return None

    def _assign_merger_type(self, m):
        if m.mtype != motif_type.HELIX:
            return MotifMergerType.SPECIFIC_SEQUENCE
        else:
            if m.name[0:5] == "HELIX":
                return MotifMergerType.NON_SPECIFIC_SEQUENCE
            else:
                return MotifMergerType.SPECIFIC_SEQUENCE

    def _link_motifs(self, m1, m2, m1_end, m2_end):
        m1_end_nodes = self._get_end_nodes(self._chain_graph.nodes, m1_end)
        m2_end_nodes = self._get_end_nodes(self._chain_graph.nodes, m2_end)

        mm_type_1 = self._assign_merger_type(m1)
        mm_type_2 = self._assign_merger_type(m2)

        if mm_type_2 == MotifMergerType.NON_SPECIFIC_SEQUENCE and \
           mm_type_1 == MotifMergerType.SPECIFIC_SEQUENCE:
            self._link_chains(m1_end_nodes, m2_end_nodes, mm_type_1, mm_type_2)
            self._bp_overrides[m2_end.uuid] = m1_end.uuid
        else:
            self._link_chains(m2_end_nodes, m1_end_nodes, mm_type_2, mm_type_1)
            self._bp_overrides[m1_end.uuid] = m2_end.uuid

    def _link_chains(self, dominant_nodes, auxiliary_nodes, mm_type_1, mm_type_2):
        if dominant_nodes[0] == dominant_nodes[1]:
            self._connect_chains(dominant_nodes[0], auxiliary_nodes[0], 1, 0,
                                 mm_type_1, mm_type_2)
            self._connect_chains(dominant_nodes[0], auxiliary_nodes[1], 0, 1,
                                 mm_type_1, mm_type_2)

        elif auxiliary_nodes[0] == auxiliary_nodes[1]:
            #print "did this happen check to make sure it worked!!!"
            self._connect_chains(dominant_nodes[1], auxiliary_nodes[0], 1, 0,
                                 mm_type_1, mm_type_2)
            self._connect_chains(dominant_nodes[0], auxiliary_nodes[0], 0, 1,
                                mm_type_1, mm_type_2)

        else:
            self._connect_chains(dominant_nodes[1], auxiliary_nodes[0], 1, 0,
                                mm_type_1, mm_type_2)
            self._connect_chains(dominant_nodes[0], auxiliary_nodes[1], 0, 1,
                                mm_type_1, mm_type_2)

    def _connect_chains(self, d_node, a_node, d_i, a_i, mm_type_1, mm_type_2):
        if a_i == 0:
            a_node.data.prime5_override = 1
            if mm_type_1 == MotifMergerType.SPECIFIC_SEQUENCE and \
               mm_type_2 == MotifMergerType.SPECIFIC_SEQUENCE:
                if a_node.data.c.first().name != d_node.data.c.last().name:
                    user_warnings.warn(
                        'overriding residues of two different types, this is likely to '
                        'produce a merged structure that is wrong!!\n',
                        user_warnings.MotifMergerWarning)

            self._res_overrides[a_node.data.c.first().uuid] = \
                d_node.data.c.last().uuid
        else:
            a_node.data.prime3_override = 1
            if mm_type_1 == MotifMergerType.SPECIFIC_SEQUENCE and \
               mm_type_2 == MotifMergerType.SPECIFIC_SEQUENCE:
                if a_node.data.c.last().name != d_node.data.c.first().name:
                    user_warnings.warn(
                        'overriding residues of two different types, this is likely to '
                        'produce a merged structure that is wrong!!\n',
                        user_warnings.MotifMergerWarning)

            self._res_overrides[a_node.data.c.last().uuid] = \
                d_node.data.c.first().uuid
        self._chain_graph.connect(d_node.index, a_node.index, d_i, a_i)

    def _get_end_nodes(self, nodes, end):
        end_nodes = [None, None]

        for n in nodes:
            for r_uuid in end.res_uuids():
                if n.data.c.first().uuid == r_uuid and end_nodes[0] is None:
                    end_nodes[0] = n
                    continue
                elif n.data.c.first().uuid == r_uuid:
                    if len(end_nodes[0].data.c) == 1 and end_nodes[1] is None:
                        end_nodes[1] = end_nodes[0]
                        end_nodes[0] = n
                        continue
                    elif len(n.data.c) == 1:
                        pass
                    else:
                        raise ValueError("cannot build chain map two residues are assigned"
                                         "to 5' chain")

                if n.data.c.last().uuid == r_uuid and end_nodes[1] is None:
                    end_nodes[1] = n
                elif n.data.c.last().uuid == r_uuid:
                    if len(end_nodes[1].data.c) == 1 and end_nodes[0] is None:
                        end_nodes[0] = end_nodes[1]
                        end_nodes[1] = n
                    else:
                        raise ValueError("cannot build chain map two residues are assigned"
                                         "to 3' chain")

        if end_nodes[0] is None or end_nodes[1] is None:
            raise ValueError("did not build map properly, both chains are not found")

        return end_nodes

    def _build_structure(self):
        raise ValueError("not implemented")

class MotifMerger(Merger):

    __slots__ = [
        "_all_bps",
        "_motifs",
        "_res_overrides",
        "_bp_overrides",
        "_rebuild_structure",
        "_chain_graph",
        "_rna_structure",
        "_aligned",
        "_m_pos",
        "_mf"
    ]

    def __init__(self, mf):
        super(self.__class__, self).__init__()
        self._mf = mf

    def __get_ordered_ends(self, s, ends):
        if len(ends) == 0:
            return []

        ordered_infos = []
        for end in ends:
            org_m = None
            for m in self._motifs.values():
                if m.get_end(bp_uuid=end.uuid):
                    org_m = m
                    break
            end_pos = m.get_end_index(end.name)

            ordered_info = { 'end' : end, 'm' : org_m,
                             'end_pos' : end_pos,
                             'm_aligned' : self._aligned[org_m.uuid],
                             'm_pos' : self._m_pos[org_m.uuid]}
            ordered_infos.append(ordered_info)

        seen = {}
        unaligned = []
        for ordered_info in ordered_infos:
            if ordered_info['m_aligned'] == 0 and \
               ordered_info['end_pos'] == ordered_info['m'].block_end_add:
                unaligned.append(ordered_info)

        unaligned.sort(key=lambda x: x['m_pos'])
        ordered_ends = []
        for info in unaligned:
            seen[info['end']] = 1
            ordered_ends.append(info['end'])

        other_infos = []
        for info in ordered_infos:
            if info['end'] in seen:
                continue
            if len(ordered_ends) > 0:
                info['dist'] = ordered_ends[0].diff(info['end'])
            else:
                info['dist'] = 0
            other_infos.append(info)


        other_infos.sort(key=lambda x: x['dist'])
        for info in other_infos:
            ordered_ends.append(info['end'])

        return ordered_ends

    def _build_structure(self):
        starts = []

        for n in self._chain_graph.nodes:
            if n.available_pos(0):
                starts.append(n)

        chains = []
        for n in starts:
            res = []
            cur = n
            while cur is not None:
                res.extend(cur.data.included_res())
                if cur.available_pos(1):
                    cur = None
                else:
                    con = cur.connections[1]
                    cur = con.partner(cur.index)

            c = chain.Chain(res)
            chains.append(c)

        res = []
        chain_cuts = []
        for i, c in enumerate(chains):
            for r in c:
                res.append(r)
            if i < len(chains) - 1:
                chain_cuts.append(len(res))
        chain_cuts.append(len(res))

        s = structure.Structure(res, chain_cuts)
        uuids = [r.uuid for r in res]

        current_bps = []
        for bp in self._all_bps.values():
            if bp.res1_uuid in uuids and bp.res2_uuid in uuids:
                current_bps.append(bp)

        ends = ends_from_basepairs(s, current_bps)
        for end in ends:
            current_bps.remove(end)

        ordered_ends = self.__get_ordered_ends(s, ends)

        self._rna_structure = self._mf.rna_structure_from_element(
                                s, current_bps, ordered_ends, "assembled")

        self._rebuild_structure = 0


class MotifStateMerger(Merger):
    def __init__(self):
        super(self.__class__, self).__init__()

    def __get_ordered_ends(self, s, ends):
        if len(ends) == 0:
            return []

        ordered_infos = []
        for end in ends:
            org_m = None
            for m in self._motifs.values():
                if m.get_end(bp_uuid=end.uuid):
                    org_m = m
                    break
            end_pos = m.get_end_index(end.name)

            ordered_info = {'end': end, 'm': org_m,
                            'end_pos': end_pos,
                            'm_aligned': self._aligned[org_m.uuid],
                            'm_pos': self._m_pos[org_m.uuid]}
            ordered_infos.append(ordered_info)

        seen = {}
        unaligned = []
        for ordered_info in ordered_infos:
            if ordered_info['m_aligned'] == 0 and \
                            ordered_info['end_pos'] == ordered_info['m'].block_end_add:
                unaligned.append(ordered_info)

        unaligned.sort(key=lambda x: x['m_pos'])
        ordered_ends = []
        for info in unaligned:
            seen[info['end']] = 1
            ordered_ends.append(info['end'])

        other_infos = []
        for info in ordered_infos:
            if info['end'] in seen:
                continue
            if len(ordered_ends) > 0:
                info['dist'] = ordered_ends[0].diff(info['end'])
            else:
                info['dist'] = 0
            other_infos.append(info)

        other_infos.sort(key=lambda x: x['dist'])
        for info in other_infos:
            ordered_ends.append(info['end'])

        return ordered_ends

    def _build_structure(self):
        starts = []

        for n in self._chain_graph.nodes:
            if n.available_pos(0):
                starts.append(n)

        chains = []
        for n in starts:
            res = []
            cur = n
            while cur is not None:
                res.extend(cur.data.included_res())
                if cur.available_pos(1):
                    cur = None
                else:
                    con = cur.connections[1]
                    cur = con.partner(cur.index)

            c = chain.Chain(res)
            chains.append(c)

        res = []
        chain_cuts = []
        for i, c in enumerate(chains):
            for r in c:
                res.append(r)
            if i < len(chains) - 1:
                chain_cuts.append(len(res))
        chain_cuts.append(len(res))

        s = motif_state.Structure(res, chain_cuts)
        uuids = [r.uuid for r in res]

        current_bps = []
        for bp in self._all_bps.values():
            if bp.res1_uuid in uuids and bp.res2_uuid in uuids:
                current_bps.append(bp)

        ends = ends_from_basepairs(s, current_bps)
        for end in ends:
            current_bps.remove(end)

        ordered_ends = self.__get_ordered_ends(s, ends)

        end_ids = []
        for end in ordered_ends:
            end_ids.append(assign_end_id(s, current_bps, ordered_ends, end,))

        seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])

        self._rna_structure = motif_state.Motif(s, current_bps, ordered_ends,
                                                end_ids, "assembled", motif_type.UNKNOWN,
                                                0, dot_bracket)
        self._rebuild_structure = 0




















