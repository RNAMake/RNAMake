import secondary_structure
import secondary_structure_graph
import graph
import motif_type
import exceptions


class SecondaryStructureChainGraph(object):
    """
    A graph data structure to store connectivity between residues based on
    their sequence and dot_bracket notation

    Nodes have 3 possible connections connection 0 is to the residue before the
    current one going from 5' to 3'. So if the current node has the first
    residue in it would have no connection at position 0. Likewise if it was the
    last residue in a chain it would have no connection at position 1. All
    other residues will have always have both connection 0 and 1 filled. Just to
    clarify all non-paired nodes can have more then one residue in them. The
    node will has as many residues that are unpaired in it.

    Where connections are located

    ^ 5' (0)\n
    \|\n
    N -- pair (2)\n
    \|\n
    v 3' (1)\n

    Example

    +----------+-+-+----+-+-+--------+-+-+-+-+
    | Sequence |G|G|AA  |C|C|UUCG    |G|G|C|C|
    +----------+-+-+----+-+-+--------+-+-+-+-+
    | Structure|(|(|. . |(|(|. . . . |)|)|)|)|
    +----------+-+-+----+-+-+--------+-+-+-+-+
    | Node     |0|1| 2  |3|4|   5    |6|7|8|9|
    +----------+-+-+----+-+-+--------+-+-+-+-+

    """

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

    def __len__(self):
        return len(self.graph)

    def add_chain(self, data, parent_index=-1, orphan=0):
        parent = self.graph.last_node
        if parent_index != -1:
            parent = self.graph.get_node(parent_index)

        if parent is None:
            return self.graph.add_data(data, n_children=3)

        return self.graph.add_data(data, parent_index,
                                   parent_pos=1, child_pos=0,
                                   n_children=3, orphan=orphan)

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
    PAIRED = 1


class NodeData(object):
    def __init__(self, residues, type):
        self.residues, self.type = residues, type


# TODO add the ability to parse "{" and "}"
class SecondaryStructureParser(object):
    """
    Parses sequence and dot bracket notation for secondary structure. Can
    can take parsed information and build many different data structures with
    it.

    :attributes:

    `structure`: secondary_structure.Structure
        secondary_structure structure to hold sequence and secondary
        structure.
    `residues`: list of secondary_structure.Residues
        residues from structure kept as seperate variable for quick referencing
    `pairs`: list of secondary_structure.Basepairs
        keeps track of all pairs found in secondary structure

    :examples:

    ..  code-block:: python

        >>> from rnamake import secondary_structure_parser
        # get parser
        >>> p = secondary_structure_parser.SecondaryStructureParser()

        # most simple parse, builds a SecondaryStructureChainGraph which
        # keeps track of the position of each element, its sequence, what
        # its connected to and what its paired too.
        >>> g = p.parse("GG+CC", "((+))")
        >>> print g
        index: 0 sequence:G 5prime: N 3prime: 1 pair: 3
        index: 1 sequence:G 5prime: 0 3prime: N pair: 2
        index: 2 sequence:C 5prime: N 3prime: 3 pair: 1
        index: 3 sequence:C 5prime: 2 3prime: N pair: 0

        # can also parse into more complex things
        >>> p.reset()
        >>> m = p.parse_to_motif("GAG+CC", "(.(+))")
        >>> print m
        <secondary_structure.Motif( GAG&CC (.(&)) )

    """

    def __init__(self):
        self.structure = secondary_structure.Structure()
        self.residues = []
        self.pairs = []

    def reset(self):
        """
        resets class so it can be used again afterwards.

        :return: None
        """

        self.structure = secondary_structure.Structure()
        self.residues = []
        self.pairs = []

    def _setup(self, sequence, dot_bracket, structure):
        """
        helper class for checking to make sure arguments to :func:`parse` are
        useable.

        :param sequence: sequence of secondary structure to parse
        :param dot_bracket: corresponding structure of given sequence
        :param structure: structure object that can be supplied instead of
            sequence and dot_bracket

        :type sequence: str
        :type dot_bracket: str
        :type structure: secondary_structure.Structure

        :return: None
        """

        # catch nothing supplied
        if sequence is None and dot_bracket is None and structure is None:
            raise exceptions.SecondaryStructureParserException(
                "no arguments supplied to parser, need to supply sequence and "
                "dot_bracket or a existing secondary_structure.Structure")

        # supplied sequence and dot_bracket
        elif sequence is not None and dot_bracket is not None:
            try:
                self.structure = secondary_structure.Structure(
                    sequence=sequence, dot_bracket=dot_bracket)

            # incorrect secondary structure
            except exceptions.SecondaryStructureException as e:
                raise exceptions.SecondaryStructureParserException(
                    "cannot parse secondary structure, error in "
                    "secondary_structure:" + e.message)

        # catch mismatch of arguments
        elif (sequence is not None and dot_bracket is None) or \
                (sequence is None and dot_bracket is not None):
            raise exceptions.SecondaryStructureParserException(
                "incorrect arguments supplied, must suppled BOTH sequence and "
                "dot_bracket secondary structure")

        else:
            self.structure = structure

        self.residues = self.structure.residues()

    def _add_unpaired_residues_to_graph(self, g, res, is_start_res):
        """
        adds a stretch of unpaired residues to the graph, a helper function
        to :func:`parser`

        :param g: the current secondary setructure graph being built
        :param res: unpaired residues to be added to graph
        :param is_start_res: whether this is a start of a chain

        :type g: SecondaryStructureChainGraph
        :type res: list of secondary_structure.Residues
        :type is_start_res: int

        :return: None
        """

        parent_index = g.get_node_by_res(self._previous_res(res[0]))
        new_data = NodeData(res, NodeType.UNPAIRED)
        g.add_chain(new_data, parent_index, is_start_res)

    def _add_paired_res_to_graph(self, g, r, is_start_res):
        """
        adds a basepairing residue to the graph, a helper function to
        :func:`parser`

        :param g: the current secondary setructure graph being built
        :param r: current paired residue
        :param is_start_res: whether this residue is the start of a chian

        :type g: SecondaryStructureChainGraph
        :type res: secondary_structure.Residue
        :type is_start_res: int

        :return: None
        """

        pair_res = self._get_bracket_pair(r)
        new_data = NodeData([r], NodeType.PAIRED)
        parent_index = g.get_node_by_res(self._previous_res(r))
        self.pairs.append(secondary_structure.Factory.get_basepair(r, pair_res))
        g.add_chain(new_data, parent_index, is_start_res)

    def _get_previous_pair(self, r):
        """
        finds the basepair created when the corresponding "(" element was
        detected.

        :param r: residue with a ")" structure
        :type r: secondary_structure.Residue

        :return: the basepair that the residue is apart of
        :rtype: secondary_structure.Basepair
        """

        pair = None
        for p in self.pairs:
            if p.res2 == r:
                pair = p
                break

        if pair is None:
            raise exceptions.SecondaryStructureParserException(
                "cannot parse secondary structure: \n%s\n%s\n"
                % (self.structure.sequence(), self.structure.dot_bracket()) + \
                "position: %d has no matching pair " % r.num)

        return pair

    def parse(self, sequence=None, dot_bracket=None, structure=None):
        """
        parses a sequence and secondary structure into a
        :class:`SecondaryStructureChainGraph`. This graph can then easily be
        converted into to other data structures.

        :param sequence: sequence of secondary structure to parse
        :param dot_bracket: corresponding structure of given sequence
        :param structure: structure object that can be supplied instead of
            sequence and dot_bracket

        :type sequence: str
        :type dot_bracket: str
        :type structure: secondary_structure.Structure

        :return: parsed structure
        :rtype: SecondaryStructureChainGraph
        """

        # checks to make sure arguments supplied are cogent.
        self._setup(sequence, dot_bracket, structure)

        g = SecondaryStructureChainGraph()

        res = []
        self.pairs = []
        for i, r in enumerate(self.residues):
            is_start_res = self._start_of_chain(r)
            if r.dot_bracket == ".":
                res.append(r)

            elif r.dot_bracket == "(":
                if len(res) > 0:
                    self._add_unpaired_residues_to_graph(g, res, is_start_res)
                    res = []

                self._add_paired_res_to_graph(g, r, is_start_res)

            elif r.dot_bracket == ")":
                if len(res) > 0:
                    self._add_unpaired_residues_to_graph(g, res, is_start_res)
                    res = []

                pair = self._get_previous_pair(r)

                new_data = NodeData([r], NodeType.PAIRED)
                parent_index = g.get_node_by_res(self._previous_res(r))
                pos = g.add_chain(new_data, parent_index, is_start_res)
                pair_res_pos = g.get_node_by_res(pair.res1)
                g.pair_res(pair_res_pos, pos)

        if len(res) > 0:
            self._add_unpaired_residues_to_graph(g, res, 0)

        return g

    def parse_to_motifs(self, sequence=None, dot_bracket=None, structure=None):
        """
        parses secondary structure into individual motifs. Helices are
        individual basepair steps of 2 basepairs. Currently will not correctly
        parse single stranded motifs on the 5' and 3' end

        :param sequence: sequence of secondary structure to parse
        :param dot_bracket: corresponding structure of given sequence
        :param structure: structure object that can be supplied instead of
            sequence and dot_bracket

        :type sequence: str
        :type dot_bracket: str
        :type structure: secondary_structure.Structure

        :return: motifs contained in secondary structure
        :rtype: list of secondary_structure.Motifs
        """

        g = self.parse(sequence, dot_bracket, structure)
        motifs = []
        self.seen_nodes = []
        for i, n in enumerate(g.graph):
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
        """
        parses secondary structure into a graph of connected secondary structure
        motif objects. Useful if one needs to perserve the connection between
        motifs.

        :param sequence: sequence of secondary structure to parse
        :param dot_bracket: corresponding structure of given sequence
        :param structure: structure object that can be supplied instead of
            sequence and dot_bracket

        :type sequence: str
        :type dot_bracket: str
        :type structure: secondary_structure.Structure

        :return: secondary structure motif graph with all motifs in it
        :rtype: :class:`secondary_structure_graph.SecondaryStructureGraph`
        """

        motifs = self.parse_to_motifs(sequence, dot_bracket, structure)
        ssg = secondary_structure_graph.SecondaryStructureGraph()

        start_m = motifs.pop(0)
        seen_motifs = {}
        open_motifs = [[start_m, -1, -1, 0]]
        while len(open_motifs) > 0:
            current, parent_index, parent_end_index, current_end_index = open_motifs.pop(0)
            seen_motifs[current] = 1

            pos = ssg.add_motif(current, parent_index, parent_end_index,
                                m_end_index=current_end_index, new_uuids=0)
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
        """
        parses secondary structure into a secondary structure motif

        :param sequence: sequence of secondary structure to parse
        :param dot_bracket: corresponding structure of given sequence
        :param structure: structure object that can be supplied instead of
            sequence and dot_bracket

        :type sequence: str
        :type dot_bracket: str
        :type structure: secondary_structure.Structure

        :return: secondary structure motif graph with all motifs in it
        :rtype: :class:`secondary_structure.Motif`
        """

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

        # print current, last_node, len(chain.residues)

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

        p = secondary_structure.Factory.get_pose(struct, bps, ends)

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

        if len(m.residues()) == 4 and len(m.basepairs) == 2:
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
        chain, next_n = self._walk_nodes(start)
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
            return self.residues[i - 1]

    def _next_res(self, r):
        i = self.residues.index(r)
        if i == len(self.residues):
            return None
        else:
            return self.residues[i + 1]

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

        raise exceptions.SecondaryStructureParserException("cannot find pair")

