import graph
import ss_tree
import copy

class MotifTreeTopology(object):
    def __init__(self, sstree):
        self.graph = graph.GraphStatic()
        self.bp_step_size = 2

        self.seen = {}
        index = 0
        mtt_node = None

        for n in sstree:
            if n.data.type == ss_tree.SS_Type.SS_BP:
                nodes = self._get_bp_nodes(n)
                if len(nodes) < self.bp_step_size:
                    continue
                mtt_data = self._build_motif_topology_node(nodes, sstree,
                                                           MotifTopologyType.BP_STEP)

            elif n.data.type != ss_tree.SS_Type.SS_SEQ_BREAK:
                nodes = [ n.parent, n ]
                for n in n.children:
                    nodes.append(n)
                mtt_data = self._build_motif_topology_node(nodes, sstree,
                                                           MotifTopologyType.NWAY)

            else:
                continue

            if len(self.graph) == 0:
                self.graph.add_data(mtt_data, -1, -1, -1, len(mtt_data.ends))
            else:
                self._add_mtt_node(mtt_data)

    def _add_mtt_node(self, mtt_data):
        for n in self.graph:
            for i, end_1 in enumerate(n.data.ends):
                for j, end_2 in enumerate(mtt_data.ends):
                    if   end_1.bounds[0] == end_2.bounds[0] and \
                         end_2.bounds[1] == end_2.bounds[1]:
                        self.graph.add_data(mtt_data, n.index, i, j, len(mtt_data.ends))
                        return

                    elif end_1.bounds[0] == end_2.bounds[0] and \
                         end_2.bounds[1] == end_2.bounds[1]:
                        self.graph.add_data(mtt_data, n.index, i, j, len(mtt_data.ends))
                        return

        raise ValueError("could not add mtt_data to tree")

    def _get_bp_nodes(self, n):
        current = n
        nodes = []
        while current.data.type == ss_tree.SS_Type.SS_BP:
            nodes.append(current)
            if current.index in self.seen:
                break
            current = current.parent
            if current is None:
                break

        nodes = nodes[::-1]
        return nodes

    def _build_motif_topology_node(self, nodes, sstree, type):
        ss_data = sstree.seq_from_nodes(nodes)

        mtt_node = MotifTopology(type)
        for n in nodes:
            self.seen[n.index] = 1
            if n.data.type != ss_tree.SS_Type.SS_BP:
                continue
            node_ss_data = [ None, None]
            for i, ss_d1 in enumerate(n.data.ss_data):
                for j, ss_d2 in enumerate(ss_data):
                    if   ss_d1.bounds[0] == ss_d2.bounds[0]:
                        node_ss_data[0] = copy.deepcopy(ss_d2)
                        correct_ss = ""
                        for e in node_ss_data[0].ss:
                            if e == ")":
                                correct_ss += "("
                            else:
                                correct_ss += e
                        node_ss_data[0].ss = correct_ss
                        break
                    elif ss_d1.bounds[0] == ss_d2.bounds[1]:
                        node_ss_data[1] = copy.deepcopy(ss_d2)
                        correct_ss = ""
                        for e in node_ss_data[1].ss:
                            if e == "(":
                                correct_ss += ")"
                            else:
                                correct_ss += e
                        node_ss_data[1].ss = correct_ss
                        break
            if len(node_ss_data) != 2:
                return ValueError("did not find the correct number of chains for basepair")

            bounds = [n.data.ss_data[0].bounds[0], n.data.ss_data[1].bounds[1]]
            end = MotifTopologyEnd(node_ss_data, bounds)
            mtt_node.add_end(end)
        return mtt_node

    def __len__(self):
        return len(self.graph)

    def __iter__(self):
        self.graph.__iter__()
        return self

    def next(self):
        return self.graph.next()

    def starting_nodes(self):
        pass


class MotifTopologyType(object):
    BP_STEP  = 0
    TWOWAY   = 1
    NWAY     = 2
    HAIRPIN  = 3
    TCONTACT = 4


class MotifTopologyEnd(object):
    def __init__(self, ss_data, bounds):
        self.ss_data, self.bounds = ss_data, bounds

    def get_id(self):
        name = ""
        for i, ss_d in enumerate(self.ss_data):
            name += ss_d.seq + "_"
            for e in ss_d.ss:
                if e == "(":
                    name += "L"
                elif e == ")":
                    name += "R"
                else:
                    name += "U"

            if i != len(self.ss_data)-1:
                name += "_"
        return name

class MotifTopology(object):
    def __init__(self, type):
        self.type, self.ends = type, []

    def add_end(self, end):
        self.ends.append(end)
