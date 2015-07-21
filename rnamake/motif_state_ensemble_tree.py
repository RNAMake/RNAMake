import tree
import motif_state_tree
import resource_manager as rm


class MotifStateEnsembleTree(object):
    def __init__(self, mt=None):
        self.tree = tree.TreeStatic()
        if mt is not None:
            self._setup_from_mt(mt)

    def add_ensemble(self, ensemble, parent_index=-1, parent_end_index=-1):
        parent = self.tree.last_node
        if parent_index != -1:
            parent = self.tree.get_node(parent_index)

        if parent is None:
            return self.tree.add_data(ensemble,
                                      len(ensemble.members[0].motif_state.end_states))

        avail_pos = self.tree.get_available_pos(parent, parent_end_index)

        for p in avail_pos:
            if p == parent.data.block_end_add:
                continue

            return self.tree.add_data(ensemble,
                                      len(ensemble.members[0].motif_state.end_states),
                                      parent.index,
                                      p)

    def to_mst(self):
        mst = motif_state_tree.MotifStateTree()

        for i, n in enumerate(self.tree):
            state = n.data.members[0].motif_state
            if i == 0:
                mst.add_state(state)
                continue

            parent_index = n.parent_index()
            parent_end_index = n.parent_end_index()

            i = mst.add_state(state, parent_index, parent_end_index)
            if i == -1:
                raise ValueError("can not build motif state tree from mset")

        return mst

    def _setup_from_mt(self, mt):
        for i, n in enumerate(mt.graph.nodes):
            mse = rm.manager.get_motif_state_ensemble(name=n.data.end_ids[0])
            if i == 0:
                self.add_ensemble(mse)
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
                    raise ValueError("did not convert motif_tree to motif state tree properly")
                self.add_ensemble(mse, parent_index, parent_end_index)


    def __len__(self):
        return len(self.tree)

    def __iter__(self):
        self.tree.__iter__()
        return self

    def next(self):
        return self.tree.next()

