import itertools
import math
import pandas as pd
import numpy as np

import tree
import motif_state_tree
import motif_ensemble
import motif_type
import resource_manager as rm
import transformations as t
import util


class MotifStateEnsembleTree(object):
    def __init__(self, mt=None, mst=None):
        self.tree = tree.TreeStatic()
        self.connections = []
        if mt is not None:
            self._setup_from_mt(mt)
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
            state = n.data.members[0].motif_state.copy()
            state.new_uuids()
            if i == 0:
                mst.add_state(state)
                continue

            parent_index = n.parent_index()
            parent_end_index = n.parent_end_index()

            j = mst.add_state(state, parent_index, parent_end_index)
            if j == -1:
                raise ValueError("can not build motif state tree from mset")

        for c in self.connections:
            mst.connections.append(c.copy())

        return mst

    def _setup_from_mt(self, mt):
        for i, n in enumerate(mt.tree.nodes):
            found_supplied = rm.manager.has_supplied_motif_ensemble(
                                            n.data.name, n.data.ends[0].name())
            if n.data.mtype == motif_type.HELIX and len(n.data.residues()) == 4:
                mse = rm.manager.get_motif_state_ensemble(name=n.data.end_ids[0])
            elif found_supplied:
                me = rm.manager.get_supplied_motif_ensemble(n.data.name,
                                                            n.data.ends[0].name())
                mse = me.get_state()
            else:
                mse = motif_ensemble.motif_state_to_motif_state_ensemble(n.data.get_state())

            if i == 0:
                self.add_ensemble(mse)
            else:
                self.add_ensemble(mse, n.parent_index(), n.parent_end_index())

        for c in mt.connections:
            self.connections.append(c.copy())

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


class MotifStateEnsembleTreeEnumerator(object):
    def __init__(self, mtst):
        self.mtst = mtst
        self.mst = self.mtst.to_mst()

        ranges = []
        for n in self.mtst:
            ranges.append(range(0, len(n.data.members)))

        self.combos = itertools.product(*ranges)
        self.last_combo = None

    def next(self):
        c = next(self.combos)
        if c is None:
            return None
        if self.last_combo is None:
            self.last_combo = c

        for i in range (0, len(c)):
            if c[i] == self.last_combo[i]:
                continue
            else:
                self.mst.replace_state(i, self.mtst.get_node(i).data.members[c[i]].motif_state)

    def record(self, fname="summary.txt"):
        mst = self.mtst.to_mst()

        ranges = []
        for n in self.mtst:
            ranges.append(range(0, len(n.data.members)))


        df = pd.DataFrame(columns="alpha,beta,gamma,dist,rot_dist".split(","))

        combos = itertools.product(*ranges)
        last_combo = None
        j = 0
        org = [0,0,0]
        I = np.eye(3)
        for c in combos:
            if last_combo == None:
                last_combo = c

            for i in range (0, len(c)):
                if c[i] == last_combo[i]:
                    continue
                else:
                    mst.replace_state(i, self.mtst.get_node(i).data.members[c[i]].motif_state)

            d = mst.last_node().data.cur_state.end_states[1].d
            r = mst.last_node().data.cur_state.end_states[1].r
            euler = t.euler_from_matrix(r)
            dist = util.distance(d, org)

            rot_dist =util.matrix_distance(I, r)


            df.loc[j] = [euler[0]*180/math.pi,euler[1]*180/math.pi,euler[2]*180/math.pi,dist, rot_dist]

            last_combo = c
            j += 1
            if j % 1000 == 0:
                print j


        df.to_csv("test.csv")



def build_motif_sub_for_motif_state_ensemble(org_m, new_motifs, mlib,
                                             extra_mse_file="test.dat"):

    f = open(extra_mse_file, "w")

    for i, end in enumerate(org_m.ends):
        all_ms = []
        all_scores = []

        for new_m in new_motifs:
            mi = mlib.get(name=new_m.name, end_name=new_m.ends[i].name())
            all_ms.append(mi)
            all_scores.append(1)

        me = motif_ensemble.MotifEnsemble()
        me.setup(org_m.end_ids[i], all_ms, all_scores)
        org_m_key = org_m.name + "-" + end.name()

        f.write(org_m_key + "!!" + me.to_str() + "\n")
    f.flush()
    f.close()