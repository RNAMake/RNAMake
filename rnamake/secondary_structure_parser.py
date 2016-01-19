import secondary_structure
import secondary_structure_graph
import graph
import motif_type

class SecondaryStructureChainGraph(object):
    def __init__(self):
        self.graph = graph.GraphStatic()

    def __repr__(self):
        s = ""
        for n in self.graph:
            seq = ""
            for r in n.data.residues:
                seq += r.name
            conn_str = "5prime: "
            if n.connections[0] is not None:
                conn_str += str(n.connections[0].partner(n.index).index) + " "
            else:
                conn_str += "N "
            conn_str += "3prime: "
            if n.connections[1] is not None:
                conn_str += str(n.connections[1].partner(n.index).index) + " "
            else:
                conn_str += "N "

            if n.connections[2] is not None:
                conn_str += "pair: " + str(n.connections[2].partner(n.index).index)

            s += "index: " + str(n.index) + " sequence:" + seq + " " + str(conn_str) + "\n"

        return s

    def add_chain(self, data, parent_index=-1, orphan=0):
        parent = self.graph.last_node
        if parent_index != -1:
            parent = self.graph.get_node(parent_index)

        if parent is None:
            return self.graph.add_data(data, -1, -1, -1, 3)

        return self.graph.add_data(data, parent_index, 1, 0, 3, orphan=orphan)

    def get_node_by_res(self, res):
        for i, n in enumerate(self.graph):
            for r in n.data.residues:
                if r == res:
                    return i
        return -1

    def pair_res(self, n_i, n_j):
        self.graph.connect(n_i, n_j, 2, 2)


class NodeType(object):
    UNPAIRED = 0
    PAIRED    = 1


class NodeData(object):
    def __init__(self, residues, type):
        self.residues, self.type = residues, type


class SecondaryStructureParser(object):

    def parse(self, sequence=None, dot_bracket=None, structure=None):
        if sequence is not None and dot_bracket is not None:
            self.structure = secondary_structure.Structure(sequence=sequence, dot_bracket=dot_bracket)
        else:
            self.structure = structure

        self.residues = self.structure.residues()

        g = SecondaryStructureChainGraph()

        res = []
        self.pairs = []
        for i, r in enumerate(self.residues):
            is_start_res = self._start_of_chain(r)
            if r.dot_bracket == ".":
                res.append(r)

            elif r.dot_bracket == "(":
                if len(res) > 0:
                    parent_index = g.get_node_by_res(self._previous_res(res[0]))
                    new_data = NodeData(res, NodeType.UNPAIRED)
                    g.add_chain(new_data, parent_index, is_start_res)
                    res = []
                pair_res = self._get_bracket_pair(r)
                new_data = NodeData([r], NodeType.PAIRED)
                parent_index = g.get_node_by_res(self._previous_res(r))
                self.pairs.append(secondary_structure.Basepair(r, pair_res))
                g.add_chain(new_data, parent_index, is_start_res)

            elif r.dot_bracket == ")":
                if len(res) > 0:
                    parent_index = g.get_node_by_res(self._previous_res(res[0]))
                    new_data = NodeData(res, NodeType.UNPAIRED)
                    g.add_chain(new_data, parent_index, is_start_res)
                    res = []
                pair = None
                for p in self.pairs:
                    if p.res2 == r:
                        pair = p
                        break
                new_data = NodeData([r], NodeType.PAIRED)
                parent_index = g.get_node_by_res(self._previous_res(r))
                pos = g.add_chain(new_data, parent_index, is_start_res)
                pair_res_pos = g.get_node_by_res(pair.res1)
                g.pair_res(pair_res_pos, pos)


        #if len(res) > 0:
        #    parent_index = g.get_node_by_res(self._previous_res(res[0]))
        #    new_data = NodeData(res, NodeType.UNPAIRED)
        #    g.add_chain(new_data, parent_index, is_start_res)

        return g

    def parse_to_motifs(self, sequence=None, dot_bracket=None, structure=None):
        g = self.parse(sequence, dot_bracket, structure)
        motifs = []
        self.seen_nodes = []
        for i, n in enumerate(g.graph):
            #single stranded motifs
            #if i == 0 or n == g.graph.last_node:
            #    chain = secondary_structure.Chain(n.data.residues)


            if n.connections[1] is None:
                continue
            if n.data.type == NodeType.UNPAIRED:
                continue
            if n.data.residues[0].dot_bracket == ")":
                continue
            m = self._generate_motif(n)
            self.seen_nodes.append(n)
            if m is None:
                continue
            motifs.append(m)
        return motifs

    def parse_to_motif_graph(self, sequence=None, dot_bracket=None, structure=None):
        motifs = self.parse_to_motifs(sequence, dot_bracket, structure)
        ssg = secondary_structure_graph.SecondaryStructureGraph()

        start_m = motifs.pop(0)
        seen_motifs = {}
        open_motifs = [[start_m, -1, -1, 0]]
        while len(open_motifs) > 0:
            current, parent_index, parent_end_index, current_end_index = open_motifs.pop(0)
            seen_motifs[current] = 1

            pos = ssg.add_motif(current, parent_index, parent_end_index,
                                m_end_index=current_end_index)
            for m in motifs:
                if m in seen_motifs:
                    continue
                for i, end1 in enumerate(current.ends):
                    if i == current_end_index:
                        continue
                    for j, end2 in enumerate(m.ends):
                        if end1 == end2:
                            open_motifs.append([m, pos, i, j])

        return ssg

    def parse_to_motif(self, sequence=None, dot_bracket=None, structure=None):
        g = self.parse(sequence, dot_bracket, structure)
        struct = self.structure
        return self._build_motif(struct)

    def parse_to_pose(self, sequence=None, dot_bracket=None, structure=None):
        motifs = self.parse_to_motifs(sequence, dot_bracket, structure)
        p = self._build_pose(self.structure)
        p.motifs = motifs
        return p

    def _walk_nodes(self, n):
        bps_count = 0
        chain = secondary_structure.Chain()
        current = n
        last_node = current
        while current is not None:
            if current in self.seen_nodes:
                return None, None

            if current.data.type == NodeType.PAIRED:
                bps_count += 1
            chain.residues.extend(current.data.residues)
            last_node = current
            if current.connections[1] is not None:
                current = current.connections[1].partner(current.index)
            else:
                break

            if bps_count == 2:
                break

        #print current, last_node, len(chain.residues)

        return chain, last_node.connections[2].partner(last_node.index)

    def _build_pose(self, struct):
        res = struct.residues()
        bps = []
        for bp in self.pairs:
            if bp.res1 in res and bp.res2 in res:
                bps.append(bp)

        chain_ends = []
        for c in struct.chains:
            chain_ends.extend([c.first(), c.last()])

        ends = []
        for bp in bps:
            if bp.res1 in chain_ends and bp.res2 in chain_ends:
                ends.append(bp)

        p = secondary_structure.Pose(struct, bps, ends)

        for end in p.ends:
            p.end_ids.append(secondary_structure.assign_end_id_new(p, end))

        return p

    def _build_motif(self, struct):
        res = struct.residues()
        bps = []
        for bp in self.pairs:
            if bp.res1 in res and bp.res2 in res:
                bps.append(bp)

        chain_ends = []
        for c in struct.chains:
            chain_ends.extend([c.first(), c.last()])

        ends = []
        for bp in bps:
            if bp.res1 in chain_ends and bp.res2 in chain_ends:
                ends.append(bp)

        m = secondary_structure.Motif(struct, bps, ends)

        for end in m.ends:
            m.end_ids.append(secondary_structure.assign_end_id_new(m, end))

        if len(m.residues()) == 4:
            m.mtype = motif_type.HELIX
            spl = m.end_ids[0].split("_")
            m.name = spl[0][0] + spl[2][1] + "=" + spl[0][1] + spl[2][0]
        elif len(m.chains()) == 2:
            m.mtype = motif_type.TWOWAY
        elif len(m.chains()) == 1:
            m.mtype = motif_type.HAIRPIN
        else:
            m.mtype = motif_type.NWAY



        return m

    def _generate_motif(self, n):
        start = n
        chain, next_n= self._walk_nodes(start)
        chains = [chain]

        if len(chain.residues) > 0 and next is None:
            struct = secondary_structure.Structure(chains)
            return self._build_motif(struct)


        while next_n != start:
            chain, next_n = self._walk_nodes(next_n)
            if next_n is None:
                return None
            chains.append(chain)

        struct = secondary_structure.Structure(chains)

        if len(struct.residues()) < 3:
            return None

        return self._build_motif(struct)

    def _previous_res(self, r):
        i = self.residues.index(r)
        if i == 0:
            return None
        else:
            return self.residues[i-1]

    def _next_res(self, r):
        i = self.residues.index(r)
        if i == len(self.residues):
            return None
        else:
            return self.residues[i+1]

    def _start_of_chain(self, r):
        for c in self.structure.chains:
            if c.first() == r:
                return 1
        return 0

    def _end_of_chain(self, r):
        for c in self.structure.chains:
            if c.last() == r:
                return 1
        return 0

    def _get_bracket_pair(self, r_start):

        bracket_count = 0
        start = 0
        for i, r in enumerate(self.residues):
            if r_start == r and not start:
                start = 1
            elif start:
                pass
            else:
                continue
            if r.dot_bracket == "(":
                bracket_count += 1
            if r.dot_bracket == ")":
                bracket_count -= 1
                if bracket_count == 0:
                    return r

        raise ValueError("cannot find pair")

