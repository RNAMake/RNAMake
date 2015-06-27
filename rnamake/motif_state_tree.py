import motif
import motif_tree
import tree
import resource_manager
import util
import settings
import basic_io


class MotifStateTree(object):
    def __init__(self, mt=None):
        self.tree = tree.TreeStatic()
        self.clash_radius = settings.CLASH_RADIUS

        if mt is not None:
            self._setup_from_mt(mt)

    def _setup_from_mt(self, mt):
        for i, n in enumerate(mt):
            if i == 0:
                self.add_state(resource_manager.manager.get_state(n.data.name))
            else:
                parent_index = 1000
                parent_end_index = -1
                for c in n.connections:
                    if c is None:
                        continue
                    if c.partner(n.index).index < parent_index:
                        parent_index =  c.partner(n.index).index
                        parent_end_index = c.end_index(c.partner(n.index).index)

                if parent_index == 1000:
                    raise ValueError("did not convert motif tree to motif state tree properly")
                j = self.add_state(resource_manager.manager.get_state(n.data.name),
                               parent_index, parent_end_index)
                if j == -1:
                    raise ValueError("could not convert motif tree to motif state tree")

    def add_state(self, state, parent_index=-1, parent_end_index=-1):
        parent = self.tree.last_node
        if parent_index != -1:
            parent = self.tree.get_node(parent_index)

        if parent is None:
            n_data = NodeData(state)
            return self.tree.add_data(n_data, len(state.end_states), -1, -1)

        avail_pos = self.tree.get_available_pos(parent, parent_end_index)

        for p in avail_pos:
            if p == 0:
                continue

            n_data = NodeData(state)
            motif.get_aligned_motif_state(parent.data.cur_state.end_states[p],
                                          n_data.cur_state,
                                          n_data.ref_state)

            if self._steric_clash(n_data):
                continue

            return self.tree.add_data(n_data, len(state.end_states), parent.index, p)

        return -1

    def to_motif_tree(self):

        mt = motif_tree.MotifTree()
        for i, n in enumerate(self.tree):
            m = resource_manager.manager.get_motif(n.data.ref_state.name)
            if i == 0:
                mt.add_motif(m)
                continue

            parent_index = n.parent_index()
            parent_end_index = n.parent_end_index()

            j = mt.add_motif(m, parent_index, parent_end_index)
            if j == -1:
                raise ValueError("cannot convert mst to mt in to_motif_tree")

        return mt

    def get_node(self, i):
        return self.tree.get_node(i)

    def replace_state(self, i, new_state):
        n = self.get_node(i)
        if len(new_state.end_states) !=  len(n.data.ref_state.end_states):
            raise ValueError("attempted to replace a state with a different number of ends")

        old_state = n.data.ref_state

        n.data.ref_state = new_state
        n.data.cur_state = new_state
        for n in tree.transverse_tree(self.tree, i):
            parent = n.parent
            pei = n.parent_end_index()

            motif.get_aligned_motif_state(parent.data.cur_state.end_states[pei],
                                          n.data.cur_state,
                                          n.data.ref_state)

    def _steric_clash(self, new_data):
        for n in self.tree.nodes[::-1]:
            for b1 in n.data.cur_state.beads:
                for b2 in new_data.cur_state.beads:
                    dist = util.distance(b1, b2)
                    if dist < self.clash_radius:
                        return 1
        return 0

    def __len__(self):
        return len(self.tree)

    def __iter__(self):
        self.tree.__iter__()
        return self

    def next(self):
        return self.tree.next()


class NodeData(object):
    def __init__(self, ref_state):
        self.ref_state = ref_state
        self.cur_state = ref_state.copy()

