import base
import option
import motif
import motif_tree
import tree
import resource_manager as rm
import util
import settings
import basic_io
import motif_tree_topology
import motif_connection
from collections import namedtuple

def motif_state_tree_from_topology(mtt, sterics=1):
    mst = MotifStateTree(sterics=sterics)
    for i, n in enumerate(mtt.tree.nodes):
        #print i, n.data.motif_name, n.data.end_ss_id, n.data.parent_end_ss_id, n.parent_index()
        if n.data.motif_name != "":
            ms = rm.manager.get_state(name=n.data.motif_name,
                                     end_id=n.data.end_ss_id)
        else:
            ms = rm.manager.get_state(end_id=n.data.end_ss_id)
        if i == 0:
            mst.add_state(ms)
        else:
            n_parent = mst.get_node(n.parent_index())
            end_id = n.data.parent_end_ss_id
            parent_end_index = n_parent.data.ref_state.end_index_with_id(end_id)
            j = mst.add_state(ms, n.parent_index(), parent_end_index=parent_end_index)
            if j == -1:
                raise ValueError("was unable to build motifstatetree from topology")

    return mst

class MotifStateTree(base.Base):
    def __init__(self, mt=None, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.tree = tree.TreeStatic()
        self.clash_radius = settings.CLASH_RADIUS
        self.connections = []

        if mt is not None:
            self._setup_from_mt(mt)

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1}

        self.options = option.Options(options)
        self.constraints = {}

    def _setup_from_mt(self, mt):
        for i, n in enumerate(mt.tree.nodes):
            ms = rm.manager.get_state(name=n.data.name, end_id=n.data.end_ids[0],
                                      end_name=n.data.ends[0].name())
            #ms.update_res_uuids(n.data.residues())

            if i == 0:
                self.add_state(ms)
            else:
                parent_index = n.parent_index()
                parent_end_index = n.parent_end_index()

                j = self.add_state(ms, parent_index, parent_end_index)
                if j == -1:
                    raise ValueError("could not convert motif tree to motif state tree")

        for c in mt.connections:
            self.connections.append(c.copy())

    def add_state(self, state, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None):
        parent = self.tree.last_node
        if parent_index != -1:
            parent = self.tree.get_node(parent_index)

        if parent is None:
            n_data = NodeData(state)
            return self.tree.add_data(n_data, len(state.end_states), -1, -1)

        if parent_end_name is not None:
            parent_end = parent.data.cur_state.get_end_state(parent_end_name)
            parent_end_index = parent.data.cur_state.end_states.index(parent_end)

        avail_pos = self.tree.get_available_pos(parent, parent_end_index)

        for p in avail_pos:
            if p == parent.data.ref_state.block_end_add:
                continue

            n_data = NodeData(state)
            motif.get_aligned_motif_state(parent.data.cur_state.end_states[p],
                                          n_data.cur_state,
                                          n_data.ref_state)

            if self.option('sterics') and self._steric_clash(n_data):
                continue

            return self.tree.add_data(n_data, len(state.end_states), parent.index, p)

        return -1

    def add_mst(self, mst,  parent_index=-1, parent_end_index=-1,
                  parent_end_name=None):
        index_dict = {}
        for i, n in enumerate(mst):
            if i == 0:
                j = self.add_state(n.data.ref_state, parent_index=parent_index,
                                   parent_end_index=parent_end_index,
                                   parent_end_name=parent_end_name)
            else:
                ind = index_dict[n.parent_index()]
                pei = n.parent_end_index()
                j = self.add_state(n.data.ref_state, parent_index=ind, parent_end_index=pei)

            index_dict[n.index] = j

    def add_connection(self, i, j, i_bp_name="", j_bp_name=""):
        node_i = self.get_node(i)
        node_j = self.get_node(j)

        node_i_indexes = []
        node_j_indexes = []
        if i_bp_name != "":
            end_i = node_i.data.cur_state.get_end_state(i_bp_name)
            ei = node_i.data.cur_state.end_states.index(end_i)
            if not node_i.available_pos(ei):
                raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                                 "using bp: " + i_bp_name + "as its not available")
            node_i_indexes.append(ei)
            name_i = i_bp_name
        else:
            node_i_indexes = node_i.available_children_pos()
            node_i_indexes.remove(0)
            name_i = node_j.data.cur_state.end_names[node_j_indexes[0]]

        if j_bp_name != "":
            end_j = node_j.data.cur_state.get_end_state(j_bp_name)
            ei = node_j.data.cur_state.end_states.index(end_j)
            if not node_j.available_pos(ei):
                raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                                 "using bp: " + j_bp_name + "as its not available")
            node_j_indexes.append(ei)
            name_j = j_bp_name

        else:
            node_j_indexes = node_j.available_children_pos()
            node_j_indexes.remove(0)
            name_j = node_j.data.cur_state.end_names[node_j_indexes[0]]

        if len(node_i_indexes) > 1 or len(node_j_indexes) > 1:
            raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                             "its unclear which ends to attach")
        if len(node_i_indexes) == 0 or len(node_j_indexes) == 0:
            raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                             " one node has no available ends")

        #self.graph.connect(i, j, node_i_indexes[0], node_j_indexes[0])

        self.connections.append(motif_connection.MotifConnection(i, j, name_i, name_j))

    def to_motif_tree(self):
        #TODO add end names to avoid palendroms
        mt = motif_tree.MotifTree(sterics=self.option('sterics'))
        for i, n in enumerate(self.tree.nodes):
            if n.data.ref_state.name != "":
                if n.data.ref_state.end_names[0] != "":
                    m = rm.manager.get_motif(name=n.data.ref_state.name,
                                             end_name = n.data.ref_state.end_names[0],
                                             end_id=n.data.ref_state.end_ids[0])
                else:
                    m = rm.manager.get_motif(name=n.data.ref_state.name,
                                             end_id=n.data.ref_state.end_ids[0])
            else:
                m = rm.manager.get_motif(end_id=n.data.ref_state.end_ids[0])

            if i == 0:
                motif.align_motif(n.data.cur_state.end_states[0],
                                  m.ends[0],
                                  m)
                mt.add_motif(m)
                continue

            parent_index = n.parent_index()
            parent_end_index = n.parent_end_index()

            j = mt.add_motif(m, parent_index, parent_end_index)
            if j == -1:
                raise ValueError("cannot convert mst to mt in to_motif_tree")

        for c in self.connections:
            mt.add_connection(c.i, c.j, c.name_i, c.name_j)

        return mt

    def to_pdb(self, name="test.pdb"):
        return self.to_motif_tree().to_pdb(name)

    def secondary_structure(self):
        mt = self.to_motif_tree()
        return mt.secondary_structure()

    def designable_secondary_structure(self):
        mt = self.to_motif_tree()
        return mt.designable_secondary_structure()

    def to_pose(self):
        return self.to_motif_tree().to_pose()

    def write_pdbs(self, name="nodes"):
        self.to_motif_tree().write_pdbs(name)

    def get_node(self, i):
        return self.tree.get_node(i)

    def replace_state(self, i, new_state):
        n = self.get_node(i)
        if len(new_state.end_states) !=  len(n.data.ref_state.end_states):
            raise ValueError("attempted to replace a state with a different number of ends")

        old_state = n.data.ref_state


        n.data.ref_state = new_state
        n.data.cur_state = new_state.copy()

        for n in tree.transverse_tree(self.tree, i):
            parent = n.parent
            if parent is None:
                continue
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

    def last_node(self):
        return self.tree.last_node

    def remove_node(self, i):
        self.tree.remove_node(index=i)

        for c in self.connections:
            if c.i == i or c.j == i:
                self.connections.remove(c)

    def remove_node_level(self, level=None):
        self.tree.remove_node_level(level)

        for c in self.connections:
            if len(self.tree) > c.i or len(self.tree) > c.j:
                self.connections.remove(c)

    def next_level(self):
        self.tree.level += 1

    def __len__(self):
        return len(self.tree)

    def __iter__(self):
        self.tree.__iter__()
        return self

    def next(self):
        return self.tree.next()

    def copy(self):
        return str_to_motif_state_tree(self.topology_to_str(), sterics=0)

    def topology_to_str(self):
        s = ""
        for n in self.tree.nodes:
            s += n.data.ref_state.name + "," + n.data.ref_state.end_ids[0] + "," + \
                 str(n.parent_index()) + "," + str(n.parent_end_index())  +  " "
        s += "|"
        for c in self.connections:
            s += c.to_str() + " "
        return s

    def get_residue(self, uuid):
        for n in self.tree:
            for r in n.data.cur_state.residues:
                if r.uuid == uuid:
                    return r
        return None

class NodeData(object):
    def __init__(self, ref_state):
        self.ref_state = ref_state
        self.cur_state = ref_state.copy()

    def get_end_state(self, name):
        return self.cur_state.get_end_state(name)

def str_to_motif_state_tree(s, sterics=1):
    spl = s.split("|")
    node_strs = spl[0].split()
    mst = MotifStateTree(sterics=sterics)
    for n_str in node_strs:
        n_spl = n_str.split(",")
        if n_spl[0] != "":
            ms = rm.manager.get_state(name=n_spl[0], end_id=n_spl[1])
        else:
            ms = rm.manager.get_state(end_id=n_spl[1])

        mst.add_state(ms, int(n_spl[2]), int(n_spl[3]))

    conn_strs = spl[1].split()
    for c_str in conn_strs:
        c_spl = c_str.split(",")
        mst.connections.append(motif_connection.MotifConnection(int(c_spl[0]), int(c_spl[1]), c_spl[2], c_spl[3]))

    return mst














