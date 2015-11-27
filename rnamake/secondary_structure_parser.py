import secondary_structure
import graph

class SecondaryStructureChainGraph(object):
    def __init__(self):
        self.graph = graph.GraphStatic()

    def __len__(self):
        return len(self.graph)

    def __repr__(self):
        s = ""
        for n in self.graph:
            seq = ""
            for c in n.data:
                seq += c.sequence()
                if c != n.data[-1]:
                    seq += "&"
            conn_str = ""
            for i, c in enumerate(n.connections):
                if c is None:
                    continue
                conn_str += str(i) + " = " + str(c.partner(n.index).index) + ", "
            s += "index: " + str(n.index) + " sequence:" + seq + " connections " + str(conn_str) + "\n"

        return s

    def get_node_by_res(self, res):
        if len(self.graph) == 0:
            return -1
        for i, n in enumerate(self.graph):
            for c in n.data:
                for r in c.residues:
                    if r == res:
                        return i
        return -1

    def add_chains(self, chains, parent_index=-1, parent_end_index=-1, orphan=0):
        parent = self.graph.last_node
        if parent_index != -1:
            parent = self.graph.get_node(parent_index)

        if parent is None:
            return self.graph.add_data(chains, -1, -1, -1, len(chains)*2)

        if parent_end_index == -1:
            avail = self.graph.check_pos_is_value(parent, 1, 0)
            pos = 1
            if not avail:
                avail = self.graph.check_pos_is_value(parent, 3, 0)
                if not avail:
                    print chains
                    print parent.data
                    print self
                    raise ValueError("trying to contect to a node that already has a default connection")
                pos = 3
        else:
            avail = self.graph.check_pos_is_value(parent, parent_end_index)
            pos = parent_end_index

        return self.graph.add_data(chains, parent_index, pos, 0, len(chains)*2, orphan=orphan)

    def connect_bp_to_parent_bp(self, n_i, n_j):
        return self.graph.connect(n_i, n_j, 2, 3)

    def connect_bp_to_parent_chain(self, n_i, n_j):
        return self.graph.connect(n_i, n_j, 2, 1)

    def node_is_bp(self, i):
        n = self.graph.get_node(i)
        if len(n.data) == 2:
            return 1
        else:
            return 0

class SecondaryStructureParser(object):
    def __init__(self):
        pass

    def _add_chain_to_graph(self, chain, g, pos=0):
        res = self._previous_res(chain.residues[0])
        parent_index = g.get_node_by_res(res)
        if res.dot_bracket == "(":
            g.add_chains([chain], parent_index, parent_end_index=1)
        else:
            g.add_chains([chain], parent_index, parent_end_index=3)



        """if pos == 0:
            g.add_chains([chain], parent_index, parent_end_index=1)
        else:
            #check for hairpin
            try:
                next_index = g.get_node_by_res(self._next_res(chain.residues[-1]))
                if next_index == parent_index:
                    print g.graph.get_node(parent_index).connections
                    g.add_chains([chain], parent_index, parent_end_index=1)

            except:
                print "still made it"
                g.add_chains([chain], parent_index, parent_end_index=3)"""


    def parse(self, sequence=None, dot_bracket=None, structure=None):
        if sequence is not None and dot_bracket is not None:
            self.structure = secondary_structure.Structure(sequence=sequence, dot_bracket=dot_bracket)
        else:
            self.structure = structure

        self.residues = self.structure.residues()

        g = SecondaryStructureChainGraph()

        chain = secondary_structure.Chain()
        for i, r in enumerate(self.residues):
            if r.dot_bracket == ".":
                chain.residues.append(r)

            elif r.dot_bracket == "(":
                #add previous unpaired region
                if len(chain.residues) > 0:
                    self._add_chain_to_graph(chain, g, pos=0)
                    chain = secondary_structure.Chain()

                pair_res = self._get_bracket_pair(r)
                chains = [ secondary_structure.Chain([r]),
                           secondary_structure.Chain([pair_res]) ]
                if self._start_of_chain(r):
                    g.add_chains(chains, -1, orphan=1)
                else:
                    parent_index = g.get_node_by_res(self._previous_res(r))
                    g.add_chains(chains, parent_index)

            elif r.dot_bracket == ")":
                #add previous unpaired region
                if len(chain.residues) > 0:
                    self._add_chain_to_graph(chain, g, pos=1)
                    chain = secondary_structure.Chain()
                if self._start_of_chain(r):
                    continue
                parent_index = g.get_node_by_res(self._previous_res(r))
                node_index = g.get_node_by_res(r)
                #print parent_index, node_index, r.name, self._previous_res(r).name, len(g)
                if g.node_is_bp(parent_index):
                    g.connect_bp_to_parent_bp(node_index, parent_index)
                else:
                    g.connect_bp_to_parent_chain(node_index, parent_index)

        return g

    def _get_basepairs(self, g):
        basepairs = []
        for n in g.graph:
            if len(n.data) == 1:
                continue
            bp = secondary_structure.Basepair(n.data[0].residues[0],
                                              n.data[1].residues[0])
            basepairs.append(bp)
        return basepairs

    def parse_to_motifs(self, sequence=None, dot_bracket=None, structure=None):
        g = self.parse(sequence, dot_bracket, structure)
        print g
        exit()
        motifs = []
        basepairs = self._get_basepairs(g)
        for i, n in enumerate(g.graph):
            if n.connections[1] is None:
                continue
            if len(n.data) == 1:
                continue

            m = self._generate_motif(n, basepairs)
            motifs.append(m)
        return motifs

    def _walk_nodes(self, n, seen, d=0):
        bps_count = 0
        chain = secondary_structure.Chain()
        current = n
        last_node = current
        while current is not None:
            if len(current.data) > 1:
                bps_count += 1

            if d == 0 or len(current.data) == 1:
                chain.residues.extend(current.data[0].residues)
            else:
                if current.data[0].residues[0] not in seen:
                    chain.residues.extend(current.data[0].residues)
                    seen[current.data[0].residues[0]] = 1
                else:
                    chain.residues.extend(current.data[1].residues)

            last_node = current
            if d == 0 or len(current.data) == 1:
                if current.connections[1] is not None:
                    current = current.connections[1].partner(current.index)
                else:
                    current = None
            else:
                if current.connections[3] is not None:
                    current = current.connections[3].partner(current.index)
                else:
                    current = None

            if bps_count == 2:
                break

        return chain, last_node

    def _generate_motif(self, n, basepairs):
        start = n
        seen = {}
        chain, next_n= self._walk_nodes(start, seen, d=0)
        chains = [chain]
        for r in chain.residues:
            seen[r] = 1

        print chain
        for r in chain.residues:
            print r.num,
        print
        print next_n.data[0].residues[0].num, next_n.data[1].residues[0].num
        while next_n != start:
            chain, next_n = self._walk_nodes(next_n, seen, d=1)
            chains.append(chain)


        struct = secondary_structure.Structure(chains)
        res = struct.residues()
        bps = []
        for bp in basepairs:
            if bp.res1 in res and bp.res2 in res:
                bps.append(bp)

        chain_ends = []
        for c in chains:
            chain_ends.extend([c.first(), c.last()])

        ends = []
        for bp in bps:
            if bp.res1 in chain_ends and bp.res2 in chain_ends:
                ends.append(bp)

        m = secondary_structure.Motif(struct, bps, ends)
        return m




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

