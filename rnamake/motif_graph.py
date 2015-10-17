import base
import option
import settings
import motif
import util
import residue
import motif_type
import motif_tree_merger
import graph
import resource_manager as rm
import structure
import rna_structure
import chain
import secondary_structure
import secondary_structure_factory as ssf


class MotifGraph(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()
        self.graph = graph.GraphStatic()
        self.clash_radius = settings.CLASH_RADIUS
        self.merger = motif_tree_merger.MotifTreeMerger()
        self.structure = MotifGraphStructure()

    def copy(self):
        mg = MotifGraph()
        new_graph = self.graph.copy()
        mg.graph = new_graph
        mg.structure = self.structure.copy(new_graph)

        return mg

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1}

        self.options = option.Options(options)
        self.constraints = {}

    def add_motif(self, m=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None):
        parent = self.graph.last_node
        if parent_index != -1:
            parent = self.graph.get_node(parent_index)

        if parent is None:
            m_copy = m.copy()
            m_copy.get_beads(m_copy.ends)
            m_copy.new_res_uuids()
            pos = self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends))
            self.structure.add(self.graph.get_node(pos))
            return pos

        if parent_end_name is not None:
            parent_end = parent.data.get_basepair(name=parent_end_name)[0]
            parent_end_index = parent.data.ends.index(parent_end)

        avail_pos = self.graph.get_availiable_pos(parent, parent_end_index)

        for p in avail_pos:
            if p == parent.data.block_end_add:
                continue

            m_added = motif.get_aligned_motif(parent.data.ends[p], m.ends[0], m)
            if self.option('sterics'):
                if self._steric_clash(m_added):
                    continue

            m_added.new_res_uuids()

            pos =  self.graph.add_data(m_added, parent.index, p, 0, len(m_added.ends))
            self.structure.add(self.graph.get_node(pos))
            return pos

        #self._update_beads(parent, new_node)

        return -1

    def remove_motif(self, pos):
        n = self.graph.get_node(pos)
        self.structure.remove(n)
        self.graph.remove_node(pos)

    def _add_motif_to_graph(self, m, parent=None, parent_end_index=None):
        if parent == None:
            m_copy = m.copy()
            m_copy.get_beads(m_copy.ends)
            m_copy.new_res_uuids()
            pos  =  self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends), orphan=1)
            self.structure.add(self.graph.get_node(pos))
            return pos

        else:
            m_added = motif.get_aligned_motif(parent.data.ends[parent_end_index],
                                              m.ends[0],
                                              m)
            m_added.new_res_uuids()
            pos = self.graph.add_data(m_added, parent.index, parent_end_index,
                                       0, len(m_added.ends))
            self.structure.add(self.graph.get_node(pos))
            return pos

    def replace_ideal_helices(self):
        size = len(self.graph)
        for i in range(0, size):
            n = self.graph.get_node(i)
            if n.data.mtype != motif_type.HELIX:
                continue
            if len(n.data.residues()) == 4:
                continue

            parent = None
            parent_end_index = None
            other = None
            other_end_index = None
            if n.connections[0] != None:
                parent = n.connections[0].partner(n.index)
                parent_end_index = n.connections[0].end_index(parent.index)
            if n.connections[1] != None:
                other = n.connections[1].partner(n.index)
                other_end_index = n.connections[1].end_index(other.index)


            name_spl = n.data.name.split(".")
            count = int(name_spl[2])
            self.remove_motif(i)
            h = rm.manager.get_motif(name="HELIX.IDEAL")
            if parent is None:
                pos = self._add_motif_to_graph(h)

            else:
                pos = self._add_motif_to_graph(h, parent, parent_end_index)

            for j in range(0, count):
                pos = self._add_motif_to_graph(h, self.graph.get_node(pos), 1)

            if other:
                self.graph.connect(pos, other.index, 1, other_end_index)
                self.structure.connect_node_chains_new(self.graph.get_node(pos),
                                                   other)

    def secondary_structure(self):
        return self.structure.secondary_structure(self.graph)

    def designable_secondary_structure(self):
        ss = self.structure.secondary_structure(self.graph)

        """if self._has_ideal_helices():
            raise ValueError("make sure there are no ideal helices, "
                             "call replace_ideal_helices")"""

        for n in self.graph.nodes:
            if n.data.name != "HELIX.IDEAL":
                continue
            for r in n.data.residues():
                r_ss = ss.get_residue(uuid=r.uuid)
                if r_ss is not None:
                    r_ss.name = "N"

        return ss

    def replace_helix_sequence(self, ss):
        print ss

        for n in self.graph.nodes:
            if n.data.mtype != motif_type.HELIX:
                continue
            bp_names = []
            for i, end in enumerate(n.data.ends):
                bp_name = ""
                for r in end.residues():
                    try:
                        bp_name += ss.get_residue(uuid=r.uuid).name
                    except:
                        #pass
                        new_uuid =  ss.residue_map[r.uuid]
                        bp_name += ss.get_residue(uuid=new_uuid).name
                bp_names.append(bp_name)

            bp_names[1] = bp_names[1][::-1]
            bp_name = "=".join(bp_names)

            m = rm.manager.get_motif(name=bp_name)

            parent = None
            parent_end_index = None
            if n.connections[0] != None:
                parent = n.connections[0].partner(n.index)
                parent_end_index = n.connections[0].end_index(parent.index)
            other = None
            other_end_index = None
            if n.connections[1] != None:
                other = n.connections[1].partner(n.index)
                other_end_index = n.connections[1].end_index(other.index)

            if parent is not None:
                m_added = motif.get_aligned_motif(parent.data.ends[parent_end_index],
                                                 m.ends[0],
                                                 m)
                n.data = m_added
            else:
                n.data = m

            if other is None:
                continue
            if other.data.mtype == motif_type.HELIX:
                continue

            m_added = motif.get_aligned_motif(n.data.ends[1],
                                              other.data.ends[0],
                                              other.data)
            other.data = m_added

    def _steric_clash(self, m):
        beads = m.beads
        for n in self.graph:
            for c1 in n.data.beads:
                for c2 in beads:
                    if c1.btype == residue.BeadType.PHOS or \
                       c2.btype == residue.BeadType.PHOS:
                        continue
                    dist = util.distance(c1.center, c2.center)
                    if dist < self.clash_radius:
                        return 1
        return 0

    def _has_ideal_helices(self):
        for n in self.graph:
            if n.data.mtype == motif_type.HELIX and len(n.data.residues()) > 4:
                return 1
        return 0

    def write_pdbs(self,name="node"):
        for n in self.graph.nodes:
            n.data.to_pdb(name+"."+str(n.index)+".pdb")

    def last_node(self):
        return self.graph.last_node

    def leafs(self):
        leaf_nodes = []
        for n in self.graph:
            f_conn = 0
            for c in n.connections:
                if c is None:
                    f_conn += 1
            if f_conn == 0:
                continue
            leaf_nodes.append(n)
        return leaf_nodes


class MotifGraphStructure(rna_structure.RNAStructure):
    def __init__(self):
        super(self.__class__, self).__init__(structure.Structure(), [])
        self.residue_map = {}
        self.basepair_map = {}
        self.all_res = {}
        self.all_bp = {}

    def connect_node_chains_new(self, n, p, conn=None):
        if conn == None:
            for c in n.connections:
                if c is None:
                    continue
                if p == c.partner(n.index):
                    conn = c
                    break

        if conn is None:
            raise ValueError("could not find connection between nodes")

        end = p.data.ends[conn.end_index(p.index)]
        n_end = n.data.ends[conn.end_index(n.index)]

        s_chain_map = rna_structure.get_chain_end_map(self.chains(), end)
        n_chain_map = rna_structure.get_chain_end_map(self.chains(), n_end)

        if   p.data.mtype == motif_type.HELIX and n.data.mtype != motif_type.HELIX:
            chains = self._merge_chains(n_chain_map, s_chain_map)
            self.basepairs.remove(end)
        else:
            chains = self._merge_chains(s_chain_map, n_chain_map)

        #update chains in structure
        for c in s_chain_map.chains() + n_chain_map.chains():
            if c in self.structure.chains:
                self.structure.chains.remove(c)
        for c in chains:
            self.structure.chains.append(c)

        #update ends
        self.ends.remove(end)

        for e in n.data.ends:
            if e != n_end:
                self.ends.append(e)

        #update basepairs
        uuids = [r.uuid for r in self.residues()]
        for bp in n.data.basepairs:
            if bp.res1.uuid in uuids and bp.res2.uuid in uuids:
                self.basepairs.append(bp)

    def connect_node_chains(self, n, p, conn=None):
        if conn == None:
            for c in n.connections:
                if c is None:
                    continue
                if p == c.partner(n.index):
                    conn = c
                    break

        if conn is None:
            raise ValueError("could not find connection between nodes")

        end = p.data.ends[conn.end_index(p.index)]
        n_end = n.data.ends[conn.end_index(n.index)]

        if end not in self.ends:
            raise ValueError("cannot find structure end")

        n_chains = [c.subchain(0) for c in n.data.chains()]
        s_chain_map = rna_structure.get_chain_end_map(self.chains(), end)
        n_chain_map = rna_structure.get_chain_end_map(n_chains, n_end)

        if   p.data.mtype == motif_type.HELIX and n.data.mtype != motif_type.HELIX:
            chains = self._merge_chains(n_chain_map, s_chain_map)
            self.basepairs.remove(end)
        else:
            chains = self._merge_chains(s_chain_map, n_chain_map)

        #update chains in structure
        for c in s_chain_map.chains():
            if c in self.structure.chains:
                self.structure.chains.remove(c)
        for c in chains:
            self.structure.chains.append(c)

        #update ends
        self.ends.remove(end)

        for e in n.data.ends:
            if e != n_end:
                self.ends.append(e)

        #update basepairs
        uuids = [r.uuid for r in self.residues()]
        for bp in n.data.basepairs:
            if bp.res1.uuid in uuids and bp.res2.uuid in uuids:
                self.basepairs.append(bp)

    def add(self, n):

        for r in n.data.residues():
            self.all_res[r.uuid] = r
        for bp in n.data.basepairs:
            self.all_bp[bp.uuid] = bp


        n_conn = 0
        for c in n.connections:
            if c is not None:
                n_conn += 1

        if len(self.structure.chains) == 0 or n_conn == 0:
            self.structure.chains.extend([c.subchain(0) for c in n.data.chains() ])
            self.basepairs.extend(n.data.basepairs)
            self.ends.extend(n.data.ends)
            return

        conn = None
        for c in n.connections:
            if c is not None:
                conn = c
                break

        p = conn.partner(n.index)
        self.connect_node_chains(n, p, conn)

    def _get_key_for_value(self, uuid):
        for k,v in self.residue_map.iteritems():
            if v == uuid:
                return k
        return None

    def _remove_res_from_chains(self, n):
        chains = self.chains()
        new_chains = []
        removed_res = []
        removed = 1
        while removed:
            removed = 0
            split = []
            for c in chains:
                for c1 in n.data.chains():
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

                    new_chains.extend(next_chains)
                    removed_res.extend(found)
                    split.append(c)
                    removed = 1

            if len(new_chains) == 0:
                break

            #un used chains to still be included
            for c in chains:
                if c not in split:
                    new_chains.append(c)

            chains = new_chains
            new_chains = []

        self.structure.chains = chains
        res = self.residues()
        for bp in self.all_bp.values():
            if bp in self.basepairs:
                continue
            if bp.res1 in res and bp.res2 in res:
                self.basepairs.append(bp)

        return removed_res

    def remove(self, n):
        removed_res = self._remove_res_from_chains(n)
        keep_bps = []
        removed_bp = []
        for bp in self.basepairs:
            if bp.res1 in removed_res and bp.res2 in removed_res:
                removed_bp.append(bp)
                continue
            keep_bps.append(bp)
        self.basepairs = keep_bps
        self.ends = rna_structure.ends_from_basepairs(self.structure, self.basepairs)

        #remove stuff from dictionaries
        for r in removed_res:
            self.all_res.pop(r.uuid, None)
            if r.uuid in self.residue_map:
                self.residue_map.pop(r.uuid)
            key = self._get_key_for_value(r.uuid)
            if key is not None:
                self.residue_map.pop(key, None)

        for bp in removed_bp:
            self.all_bp.pop(bp, None)

    def _merge_chains(self, cm1, cm2):
        merged_chain_1, merged_chain_2 = None, None
        if   cm1.is_hairpin() and cm2.is_hairpin():
            raise ValueError("cannot merge an hairpin with another hairpin")
        elif cm1.is_hairpin():
            pass
        elif cm2.is_hairpin():
            pass
        else:
            merged_chain_1 = self._get_merged_chain(cm1.p5_chain, cm2.p3_chain, 1, 1)
            merged_chain_2 = self._get_merged_chain(cm1.p3_chain, cm2.p5_chain, 0, 1)

        return merged_chain_1, merged_chain_2

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

    def secondary_structure(self, graph):
        ss = ssf.factory.secondary_structure_from_motif(self)

        ss_motifs = []
        for n in graph.nodes:
            ss_chains = []
            for c in n.data.chains():
                ss_res = []
                for r in c.residues:
                    correct_r = self.get_residue(uuid=r.uuid)
                    ss_r = ss.get_residue(uuid=correct_r.uuid)
                    ss_res.append(ss_r)
                ss_chains.append(secondary_structure.Chain(ss_res))
            ss_struct = secondary_structure.Structure(ss_chains)
            ss_bps = []
            for bp in self.basepairs:
                res1 = ss_struct.get_residue(uuid=bp.res1.uuid)
                res2 = ss_struct.get_residue(uuid=bp.res2.uuid)
                if res1 is None or res2 is None:
                    continue
                ss_bp = secondary_structure.Basepair(res1, res2, bp.uuid)
                ss_bps.append(ss_bp)
            ss_rna_struct = secondary_structure.RNAStructure(
                                ss_struct, ss_bps, [], n.data.name, n.data.path,
                                n.data.mtype, n.data.score, n.data.end_ids)

            ss_ends = []
            for end in n.data.ends:
                res1 = self.get_residue(uuid=end.res1.uuid)
                res2 = self.get_residue(uuid=end.res2.uuid)
                ss_res1 = ss_struct.get_residue(uuid=res1.uuid)
                ss_res2 = ss_struct.get_residue(uuid=res2.uuid)
                ss_bp = ss_rna_struct.get_basepair(ss_res1, ss_res2)
                ss_ends.append(ss_bp)
            ss_rna_struct.ends = ss_ends
            ss_motifs.append(secondary_structure.Motif(r_struct=ss_rna_struct,
                                                       id=n.data.id))

        ss_struct = secondary_structure.Structure(ss.chains)
        ss_p = secondary_structure.Pose(ss_struct, ss.basepairs, ss.ends)
        ss_p.residue_map = self.residue_map
        ss_p.motifs = ss_motifs

        return ss_p

    def copy(self, graph):
        new_all_res = {}
        new_all_bp = {}
        for n in graph.nodes:
            for r in n.data.residues():
                new_all_res[r.uuid] = r
            for bp in n.data.basepairs:
                new_all_bp[bp.uuid] = bp


        new_struct = MotifGraphStructure()
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

        return new_struct











