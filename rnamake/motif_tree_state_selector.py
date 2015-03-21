import motif_tree_state
import motif_type
import settings
import motif_tree_precomputer
import time

class MotifTreeStateSelector(object):
    def __init__(self, mts_libs=None, mode="helix_flank", nodes=None):
        self.clash_lists = {}
        if nodes is not None:
            self.nodes = nodes
            return

        if mts_libs is None:
            raise ValueError("must specify mts_libs or nodes")

        if mode == "all":
            self.nodes = []
            for i, mts_lib in enumerate(mts_libs):
                self.nodes.append(SelectorNode(mts_lib, i))
            for n1 in self.nodes:
                for n2 in self.nodes:
                    n1.add_connection(n2)

        elif mode == "helix_flank":
            hlib = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
            self.nodes = [SelectorNode(hlib, 0)]
            for i, mts_lib in enumerate(mts_libs):
                self.nodes.append(SelectorNode(mts_lib, i+1))
            for i, n in enumerate(self.nodes):
                if i == 0:
                    continue
                n.add_connection(self.nodes[0])
                self.nodes[0].add_connection(n)

        self._setup_clash_lists()

    def _setup_clash_lists(self):
        seen = {}
        for n in self.nodes:
            for n2 in n.connections:
                if n.mts_lib.mtype == motif_type.UNKNOWN:
                    continue
                if n2.mts_lib.mtype == motif_type.UNKNOWN:
                    continue
                name = settings.PRECOMPUTED_PATH + "motif_tree_states/"
                name += motif_type.type_to_str(n.mts_lib.mtype) + "_"
                name += motif_type.type_to_str(n2.mts_lib.mtype) + ".clist"
                if name in seen:
                    continue
                seen[name] = 1
                try:
                    f = open(name)
                    lines = f.readlines()
                    f.close()
                    clist = {}
                    for l in lines:
                        clist[l.rstrip()] = 2
                    self.clash_lists[str(n.mts_lib.mtype) + "-" + \
                                     str(n2.mts_lib.mtype)] = clist
                except:
                    print "warning no clist file loaded for ",
                    motif_type.type_to_str(n.mts_lib.mtype), " and ",
                    motif_type.type_to_str(n2.mts_lib.mtype)

    def get_children_mts(self, node):
        lib_type = node.lib_type
        if lib_type != -1:
            connections = self.nodes[lib_type].connections
        else:
            connections = [self.nodes[0]]
        children, types = [], []
        for c in connections:
            if c.max_uses <= node.lib_type_usage(c.index):
                continue
            clist_name = str(node.lib_type) + "-" + str(c.index)
            clist = {}
            if clist_name in self.clash_lists:
                clist = self.clash_lists[clist_name]
            for mts in c.mts_lib.motif_tree_states:
                #key = node.mts.name + " " + mts.name
                #if key in clist:
                #    continue
                children.append(mts)
                types.append(c.index)

        return children, types

    def is_valid_solution(self, current):
        for i, n in enumerate(self.nodes):
            if n.required_uses > current.node_counts[i]:
                return 0
        return 1

    def set_node_uses(self, index, max_uses=None, required_uses=None):
        if max_uses is not None:
            self.nodes[index].max_uses = max_uses
        if required_uses is not None:
            self.nodes[index].required_uses = required_uses


class SelectorNode(object):
    def __init__(self, mts_lib, index, max_uses=1000, required_uses=0):
        self.mts_lib = mts_lib
        self.index = index
        self.connections = []
        self.max_uses = max_uses
        self.required_uses = required_uses

    def add_connection(self, node):
        self.connections.append(node)


def default_selector(types=None, mode=None):
    if mode is None:
        mode = "helix_flank"

    if types is None:
        types = [motif_type.TWOWAY]
    mts_libs = []
    for t in types:
        mts_lib = motif_tree_state.MotifTreeStateLibrary(t)
        mts_libs.append(mts_lib)

    return MotifTreeStateSelector(mts_libs, mode=mode)


def force_include_motif_selector(m, types=None, mode=None):
    if mode is None:
        mode = "helix_flank"

    if types is None:
        types = [motif_type.TWOWAY]

    mtp = motif_tree_precomputer.MotifTreePrecomputer(name="temp",
                                                      max_bps_per_end=0)
    mtp.precompute_motif(m)
    mts_libs = []
    for t in types:
        mts_lib = motif_tree_state.MotifTreeStateLibrary(t)
        mts_libs.append(mts_lib)
    m_mts_lib = motif_tree_state.MotifTreeStateLibrary(libpath="temp.new.me")
    mts_libs.append(m_mts_lib)
    selector = MotifTreeStateSelector(mts_libs, mode=mode)
    selector.set_node_uses(2, 1, 1)
    return selector











