import motif_type
import rna_structure
import chain
import secondary_structure_factory as ssf
import secondary_structure
import util
import graph
import user_warnings


class ChainNodeData(object):
    def __init__(self, c, m_id, prime5_override=0, prime3_override=0):
        self.c = c
        self.m_id = m_id
        self.prime5_override = prime5_override
        self.prime3_override = prime3_override

    def included_res(self):
        res = self.c.residues[::]
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


class MotifMerger(object):
    def __init__(self):
        #super(self.__class__, self).__init__()
        self.all_bps = {}
        self.motifs = {}
        self.res_overrides = {}
        self.bp_overrides = {}

        self.rebuild_structure = 1
        self.chain_graph = graph.GraphStatic()

        self.rna_structure = rna_structure.RNAStructure()

    def get_structure(self):
        if self.rebuild_structure == 1:
            self._build_structure()
            self.rebuild_structure = 0

        return self.rna_structure

    def _build_structure(self):
        starts = []

        for n in self.chain_graph.nodes:
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

        self.rna_structure.structure.chains = chains
        res = self.rna_structure.residues()
        uuids = [r.uuid for r in res]

        current_bps = []
        for bp in self.all_bps.values():
            if bp.res1.uuid in uuids and bp.res2.uuid in uuids:
                current_bps.append(bp)

        self.rna_structure.basepairs = current_bps
        self.rna_structure.ends = rna_structure.ends_from_basepairs(
            self.rna_structure.structure, current_bps)

        self.rebuild_structure = 0

    def copy(self, new_motifs=None):
        new_merger = MotifMerger()
        new_merger.chain_graph = self.chain_graph.copy()
        new_merger.res_overrides = { k : v for k,v in self.res_overrides.iteritems() }
        new_merger.bp_overrides = {k : v for k,v in self.bp_overrides.iteritems() }

        for m in new_motifs:
            new_merger.update_motif(m)

        return new_merger

    def add_motif(self, m, m_end=None, parent=None, parent_end=None):
        new_chains = [c.subchain(0) for c in m.chains()]

        for c in new_chains:
            data = ChainNodeData(c, m.id)
            self.chain_graph.add_data(data, n_children=2, orphan=1)

        for bp in m.basepairs:
            self.all_bps[bp.uuid] = bp

        if parent is not None:
            self._link_motifs(parent, m, parent_end, m_end)

        self.motifs[m.id] = m
        self.rebuild_structure = 1

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

    def secondary_structure(self):
        ss = ssf.factory.secondary_structure_from_motif(self.get_structure())

        ss_motifs = []
        for m in self.motifs.values():
            ss_chains = []
            for c in m.chains():
                ss_res = []
                for r in c.residues:
                    r_cur = r
                    if r.uuid in self.res_overrides:
                        r_cur = self.get_residue(self.res_overrides[r.uuid])
                    ss_r = ss.get_residue(uuid=r_cur.uuid)
                    if ss_r is None:
                        raise ValueError("could not find residue during ss build: " +
                                         str(r_cur) + "from " + m.name)
                    ss_res.append(ss_r)
                ss_chains.append(secondary_structure.Chain(ss_res))
            ss_struct = secondary_structure.Structure(ss_chains)
            ss_bps = []
            for bp in m.basepairs:
                if bp.bp_type != "cW-W":
                    continue
                if not util.wc_bp(bp) and not util.gu_bp(bp):
                    continue
                correct_bp = bp
                if bp.uuid in self.bp_overrides:
                    correct_bp = self.get_basepair(self.bp_overrides[bp.uuid])
                ss_bp = ss.get_basepair(uuid=correct_bp.uuid)
                if ss_bp is None:
                    raise ValueError("could not find basepair during ss build")
                ss_bps.append(ss_bp)
            ss_rna_struct = secondary_structure.RNAStructure(
                                ss_struct, ss_bps, [], m.name, m.path,
                                m.score, m.end_ids[:])
            ss_ends = []
            for end in m.ends:
                correct_bp = end
                if end.uuid in self.bp_overrides:
                    correct_bp = self.get_basepair(self.bp_overrides[end.uuid])
                ss_bp = ss.get_basepair(uuid=correct_bp.uuid)
                if ss_bp is None:
                    print m, end
                    raise ValueError("coult not find end during ss build")
                ss_ends.append(ss_bp)
            ss_rna_struct.ends = ss_ends
            ss_motifs.append(secondary_structure.Motif(r_struct=ss_rna_struct,
                                                       id=m.id,
                                                       mtype=m.mtype))

        ss_struct = secondary_structure.Structure(ss.chains())
        ss_p = secondary_structure.Pose(ss_struct, ss.basepairs, ss.ends)
        ss_p.motifs = ss_motifs
        return ss_p

    def get_residue(self, uuid):
        for m in self.motifs.values():
            r = m.get_residue(uuid=uuid)
            if r is not None:
                return r
        return None

    def get_basepair(self, uuid):
        if uuid in self.all_bps:
            return self.all_bps[uuid]
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
        m1_end_nodes = self._get_end_nodes(self.chain_graph.nodes, m1_end)
        m2_end_nodes = self._get_end_nodes(self.chain_graph.nodes, m2_end)

        mm_type_1 = self._assign_merger_type(m1)
        mm_type_2 = self._assign_merger_type(m2)

        if mm_type_2 == MotifMergerType.NON_SPECIFIC_SEQUENCE and \
           mm_type_1 == MotifMergerType.SPECIFIC_SEQUENCE:
            self._link_chains(m1_end_nodes, m2_end_nodes, mm_type_1, mm_type_2)
            self.bp_overrides[m2_end.uuid] = m1_end.uuid
        else:
            self._link_chains(m2_end_nodes, m1_end_nodes, mm_type_2, mm_type_1)
            self.bp_overrides[m1_end.uuid] = m2_end.uuid

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
                        'produce a merged structure that is wrong!!',
                        user_warnings.MotifMergerWarning)

            self.res_overrides[a_node.data.c.first().uuid] = \
                d_node.data.c.last().uuid
        else:
            a_node.data.prime3_override = 1
            if mm_type_1 == MotifMergerType.SPECIFIC_SEQUENCE and \
               mm_type_2 == MotifMergerType.SPECIFIC_SEQUENCE:
                if a_node.data.c.last().name != d_node.data.c.first().name:
                    user_warnings.warn(
                        'overriding residues of two different types, this is likely to '
                        'produce a merged structure that is wrong!!',
                        user_warnings.MotifMergerWarning)

            self.res_overrides[a_node.data.c.last().uuid] = \
                d_node.data.c.first().uuid
        self.chain_graph.connect(d_node.index, a_node.index, d_i, a_i)

    def _get_end_nodes(self, nodes, end):
        end_nodes = [None, None]

        for n in nodes:
            for r in end.residues():
                if n.data.c.first().uuid == r.uuid and end_nodes[0] is None:
                    end_nodes[0] = n
                    continue
                elif n.data.c.first().uuid == r.uuid:
                    if len(end_nodes[0].data.c) == 1 and end_nodes[1] is None:
                        end_nodes[1] = end_nodes[0]
                        end_nodes[0] = n
                        continue
                    elif len(n.data.c) == 1:
                        pass
                    else:
                        raise ValueError("cannot build chain map two residues are assigned"
                                         "to 5' chain")

                if n.data.c.last().uuid == r.uuid and end_nodes[1] is None:
                    end_nodes[1] = n
                elif n.data.c.last().uuid == r.uuid:
                    if len(end_nodes[1].data.c) == 1 and end_nodes[0] is None:
                        end_nodes[0] = end_nodes[1]
                        end_nodes[1] = n
                    else:
                        raise ValueError("cannot build chain map two residues are assigned"
                                         "to 3' chain")

        if end_nodes[0] is None or end_nodes[1] is None:
            raise ValueError("did not build map properly, both chains are not found")

        return end_nodes
