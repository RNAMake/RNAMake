import motif_type
import settings
import resource_manager
import graph

class MotifStateSelector(object):
    def __init__(self):
        self.graph = graph.GraphDynamic()

    def add(self, lib=None, mse=None, m=None, max_uses=1000, required_uses=0):
        if lib is not None:
            ms_lib = resource_manager.manager.ms_libs[lib]
            ms_lib.load_all()
            motif_states = ms_lib.all()
            d = MotifStateSelectorNodeData(ms_lib.name, motif_states,
                                           max_uses, required_uses)
            self.graph.add_data(d)

        if m is not None:
            motif_states = [m.get_state()]
            d = MotifStateSelectorNodeData(m.name, motif_states,
                                           max_uses, required_uses)
            self.graph.add_data(d)

    def connect(self, name_i, name_j):
        i, j = -1, -1
        for n in self.graph:
            if n.data.name == name_i:
                i = n.index
            if n.data.name == name_j:
                j = n.index

        if i == -1 or j == -1:
            raise ValueError("could not find a node with that name")

        self.graph.connect(i, j)

    def get_children_ms(self, node):

        connected_nodes = []
        if node.ntype != -1:
            connections = self.graph.get_node(node.ntype).connections
            for c in connections:
                connected_nodes.append(c.partner(node.ntype))

        else:
            connected_nodes = [ self.graph.get_node(0) ]
        children, types = [], []
        for c in connected_nodes:
            if c.data.max_uses <= node.node_type_usage(c.index):
                continue
            for ms in c.data.motif_states:
                children.append(ms)
                types.append(c.index)

        return children, types

    def is_valid_solution(self, current):
        for i, n in enumerate(self.graph):
            if n.data.required_uses > current.node_type_usages[i]:
                return 0
        return 1

    def score(self, current):
        diff = 0
        for i, n in enumerate(self.graph):
            if n.data.required_uses > current.node_type_usages[i]:
                diff += (n.data.required_uses - current.node_type_usages[i])
        return diff

class MSS_RoundRobin(MotifStateSelector):
    def __init__(self):
        self.graph = graph.GraphDynamic()

    def add(self, lib=None, mse=None, max_uses=1000, required_uses=0):
        super(self.__class__, self).add(lib, mse, max_uses, required_uses)
        i = len(self.graph)-1
        print len(self.graph)
        for n in self.graph.nodes:
            self.graph.connect(n.index, i)


class MSS_HelixFlank(MotifStateSelector):
    def __init__(self):
        self.graph = graph.GraphDynamic()
        self.add("ideal_helices")

    def add(self, lib=None, mse=None, m=None, max_uses=1000, required_uses=0):
        super(self.__class__, self).add(lib, mse, m, max_uses, required_uses)
        i = len(self.graph)-1
        if i != 0:
            self.graph.connect(0, i)


class MotifStateSelectorNodeData(object):
    def __init__(self, name, motif_states, max_uses, required_uses):
        self.name, self.motif_states, self.max_uses = name, motif_states, max_uses
        self.required_uses = required_uses


def default_selector():
    selector = MSS_HelixFlank()
    selector.add('unique_twoway')
    return selector










