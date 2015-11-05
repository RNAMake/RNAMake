import tree
import motif_state_tree
import motif_ensemble
import motif_type
import resource_manager as rm


class MotifStateEnsembleTree(object):
    def __init__(self, mt=None, mst=None, extra_mes=None):
        self.tree = tree.TreeStatic()
        if mt is not None:
            self._setup_from_mt2(mt, extra_mes)
        if mst is not None:
            self._setup_from_mst(mst)

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
        return -1

    def to_mst(self):
        mst = motif_state_tree.MotifStateTree(sterics=0)

        for i, n in enumerate(self.tree.nodes):
            state = n.data.members[0].motif_state
            if i == 0:
                mst.add_state(state)
                continue

            parent_index = n.parent_index()
            parent_end_index = n.parent_end_index()

            j = mst.add_state(state, parent_index, parent_end_index)
            if j == -1:
                raise ValueError("can not build motif state tree from mset")

        return mst

    def _setup_from_mt2(self, mt, extra_mes):
        if extra_mes is None:
            extra_mes = {}

        for i, n in enumerate(mt.tree.nodes):
            if n.data.mtype == motif_type.HELIX:
                mse = rm.manager.get_motif_state_ensemble(name=n.data.end_ids[0])
            else:
                if n.index in extra_mes:
                    mse = extra_mes[n.index]
                else:
                    m = n.data
                #m.name = "unknown."+str(i)+".pdb"
                #rm.manager.add_motif(motif=m)
                    mse = motif_ensemble.motif_state_to_motif_state_ensemble(m.get_state())

            mse.update_res_uuids(n.data.residues())

            if i == 0:
                self.add_ensemble(mse)
            else:
                self.add_ensemble(mse, n.parent_index(), n.parent_end_index())


    def _setup_from_mt(self, mt):
        for i, n in enumerate(mt.graph.nodes):
            try:
                mse = rm.manager.get_motif_state_ensemble(name=n.data.end_ids[0])
            except ValueError:
                m = rm.manager.get_motif(name=n.data.name, end_id=n.data.end_ids[0])
                mse = motif_ensemble.motif_state_to_motif_state_ensemble(m.get_state())

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

    def _setup_from_mst(self, mst):
        for i, n in enumerate(mst.tree.nodes):
            try:
                mse = rm.manager.get_motif_state_ensemble(name=n.data.ref_state.end_ids[0])
            except ValueError:
                mse = motif_ensemble.motif_state_to_motif_state_ensemble(n.data.ref_state)

            if i ==0:
                self.add_ensemble(mse)
            else:
                self.add_ensemble(mse, n.parent_index(), n.parent_end_index())

    def __len__(self):
        return len(self.tree)

    def __iter__(self):
        self.tree.__iter__()
        return self

    def next(self):
        return self.tree.next()

    def get_node(self, i):
        return self.tree.get_node(i)

