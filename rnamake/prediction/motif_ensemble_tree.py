import rnamake.motif_tree_state as motif_tree_state
import motif_ensemble

class MotifEnsembleTree(object):
    def __init__(self, ensemble=None):
        if ensemble is None:
            head = self._default_head()
        else:
            head = MotifEnsembleTreeNode(ensemble, None, 0)
        self.nodes = [ head ]
        self.last_node = head

    def add_ensemble(self, ensemble, parent=None, parent_end_index=None):
        if parent is None:
            parent = self.last_node
        if parent_end_index is None:
            parent_end_index = parent.available_ends()[0]
        node = MotifEnsembleTreeNode(ensemble, parent, len(self.nodes))
        parent.add_child(node, parent_end_index)
        self.nodes.append(node)
        self.last_node = node
        return node

    def _default_head(self):
        me = motif_ensemble.MotifEnsemble()
        mts = motif_tree_state.ref_mts()
        ms = motif_ensemble.MotifState(mts, 1.0)
        me.motif_states.append(ms)
        return MotifEnsembleTreeNode(me, None, 0)

    def get_mtst(self):
        mtst = None
        mts_name = self.nodes[0].motif_ensemble.motif_states[0].mts.name
        if mts_name == "start":
            mtst = motif_tree_state.MotifTreeStateTree()
        else:
            mts =  self.nodes[0].motif_ensemble.motif_states[0].mts
            mtst = motif_tree_state.MotifTreeStateTree(mts)
        open_nodes = [ self.nodes[0] ]
        while len(open_nodes) != 0:
            current = open_nodes.pop(0)
            for i, n in current.connected_nodes.iteritems():
                if n is None:
                    continue
                mts = n.motif_ensemble.motif_states[0].mts
                parent = mtst.nodes [ current.index ]
                index = parent.mts.end_indexes.index(i)
                parent_end = parent.states[ index ]
                mstn = mtst.add_state(mts, parent, parent_end)
                if mstn is None:
                    raise ValueError("cannot build mstn based on motif ensemble \
                                     tree topology")
                open_nodes.append(n)
        return mtst


class MotifEnsembleTreeNode(object):
    def __init__(self, motif_ensemble, parent, index):
        self.motif_ensemble, self.parent = motif_ensemble, parent
        self.index = index
        self.ends = self.motif_ensemble.motif_states[0].mts.end_indexes
        self.connected_nodes = {e : None for e in self.ends}

    def available_ends(self):
        ends = []
        for k, v in self.connected_nodes.iteritems():
            if v is None:
                ends.append(k)
        return ends

    def add_child(self, node, end_direction):
        if self.connected_nodes[end_direction] is not None:
            raise ValueError("cannot add child already is using this end")
        self.connected_nodes[end_direction] = node
