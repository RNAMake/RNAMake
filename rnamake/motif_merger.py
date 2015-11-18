import motif_type
import rna_structure
import chain
import secondary_structure_factory as ssf
import secondary_structure
import util
import graph

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


class MotifMerger(object):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.all_bps = {}
        self.motifs = {}

        self.rebuild_structure = 1
        self.chain_graph = graph.GraphStatic()

        self.rna_structure = rna_structure.RNAStructure()

    def get_structure(self):
        if self.rebuild_structure == 1:
            self._build_structure()
            self.rebuild_structure = 0

        return self.rna_structure

    def _build_structure(self):
        res = []
        starts = []

        for n in self.chain_graph.nodes:
            if n.available_pos(0):
                starts.append(n)

        i = 0
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
        self.rna_structure.ends = rna_structure.ends_from_basepairs(self.rna_structure.structure, current_bps)

        self.rebuild_structure = 0

    def copy(self, new_motifs=None):
        new_merger = MotifMerger()
        new_merger.chain_graph = self.chain_graph.copy()

        for m in new_motifs:
            self.update_motif(m)

        return new_merger

    def copy_old(self, new_motifs=None):
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
        new_chains = [c.subchain(0) for c in m.chains() ]

        for c in new_chains:
            data = ChainNodeData(c, m.id)
            self.chain_graph.add_data(data, n_children=2, orphan=1)

        for bp in m.basepairs:
            self.all_bps[bp.uuid] = bp

        if parent is not None:
            s_end_nodes = self._get_end_nodes(self.chain_graph.nodes, parent_end)
            m_end_nodes = self._get_end_nodes(self.chain_graph.nodes, m_end)

            self._link_chains(s_end_nodes, m_end_nodes)

        self.motifs[m.id] = m
        self.rebuild_structure = 1

    def connect_motifs(self, m1, m2, m1_end, m2_end):
        self._merge_motifs(m1, m2, m1_end, m2_end, self.chains(), self.chains())

        uuids = [r.uuid for r in self.residues()]
        for bp in m2.basepairs:
            if bp.res1.uuid in uuids and bp.res2.uuid in uuids and \
               bp not in self.basepairs:
                self.basepairs.append(bp)

    def remove_motif(self, m):
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
                if p_i == 0:
                    p.data.prime5_override = 0
                else:
                    p.data.prime3_override = 0
            self.chain_graph.remove_node(r.index)

        self.rebuild_structure = 1

    def update_motif(self, m):
        for n in self.chain_graph.nodes:
            if n.data.m_id != m.id:
                continue
            res = n.data.c.residues
            new_res = []
            for r in res:
                new_r = m.get_residue(uuid=r.uuid)
                if new_r is None:
                    raise ValueError("could not find corresponding res by uuid")
                new_res.append(new_r)
            n.data.c.residues = new_res

        for bp in m.basepairs:
            if bp.uuid in self.all_bps:
                self.all_bps[bp.uuid] = bp
        self.motifs[m.id] = m


    def update_motif_old(self, m):
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

    def _link_chains(self, dominant_nodes, auxiliary_nodes):
        if dominant_nodes[0] == dominant_nodes[1]:
            pass

        elif auxiliary_nodes[0] == auxiliary_nodes[1]:
            pass

        else:
            auxiliary_nodes[0].data.prime5_override = 1
            auxiliary_nodes[1].data.prime3_override = 1
            self.chain_graph.connect(dominant_nodes[1].index,
                                     auxiliary_nodes[0].index, 1, 0)
            self.chain_graph.connect(auxiliary_nodes[1].index,
                                     dominant_nodes[0].index, 1, 0)

            #self.links.extend([l1, l2])

    def _get_end_nodes(self, nodes, end):
        end_nodes = [None, None]

        for n in nodes:
            #for r in c.residues:
            #    print r,
            #print
            for r in end.residues():
                if n.data.c.first().uuid == r.uuid and end_nodes[0] == None:
                    end_nodes[0] = n
                elif n.data.c.first().uuid == r.uuid:
                    raise ValueError("cannot build chain map two residues are assigned"
                                     "to 5' chain")

                if n.data.c.last().uuid == r.uuid and end_nodes[1] == None:
                    end_nodes[1] = n
                elif n.data.c.last().uuid == r.uuid:
                    raise ValueError("cannot build chain map two residues are assigned"
                                     "to 3' chain")

        if end_nodes[0] == None or end_nodes[1] == None:
            raise ValueError("did not build map properly, both chains are not found")

        return end_nodes

    def _merge_chains(self, cm1, cm2):
        merged_chain_1, merged_chain_2 = None, None
        if   cm1.is_hairpin() and cm2.is_hairpin():
            raise ValueError("cannot merge an hairpin with another hairpin")
        elif cm1.is_hairpin():
            res1 = cm2.p3_chain.residues.pop()
            self.residue_map[res1.uuid] = cm1.p5_chain.residues[0].uuid
            res2 = cm2.p5_chain.residues.pop(0)
            self.residue_map[res2.uuid] = cm1.p5_chain.residues[-1].uuid
            merged_chain_1 = self._get_merged_hairpin(cm2.p3_chain, cm2.p5_chain,
                                                      cm1.p5_chain, 0, 0)
        elif cm2.is_hairpin():
            res1 = cm1.p5_chain.residues.pop(0)
            res2 = cm1.p3_chain.residues.pop()
            self.residue_map[res1.uuid] = cm2.p5_chain.residues[-1].uuid
            self.residue_map[res2.uuid] = cm2.p5_chain.residues[0].uuid

            merged_chain_1 = self._get_merged_hairpin(cm1.p3_chain, cm1.p5_chain,
                                                      cm2.p5_chain, 0, 0)
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
