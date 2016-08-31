from collections import namedtuple

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
import exceptions

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

    #SETUP FUNCTIONS ##########################################################
    def __init__(self, mt=None, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.tree = tree.TreeStatic()
        self.clash_radius = settings.CLASH_RADIUS
        self.connections = []

        if mt is not None:
            self._setup_from_mt(mt)

    def __len__(self):
        return len(self.tree)

    def __iter__(self):
        self.tree.__iter__()
        return self

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1}

        self.options = option.Options(options)
        self.constraints = {}

    def _setup_from_mt(self, mt):
        self.option('sterics', mt.option('sterics'))
        for i, n in enumerate(mt.tree.nodes):
            ms = rm.manager.get_state(name=n.data.name,
                                      end_id=n.data.end_ids[0],
                                      end_name=n.data.ends[0].name())
            ms.uuid = n.data.id
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

    def copy(self):
        return str_to_motif_state_tree(self.topology_to_str(), sterics=0)

    #ADD FUNCTIONS      #######################################################
    def _validate_arguments_to_add_state(self, ms, m_name):
        """
        makes sure the add_motif function is called correctly

        :param m: motif to add to tree
        :type m: Motif object

        :param m_name: name of motif to add
        :type m_name: str

        :return: None
        """

        if ms is not None and m_name is not None:
            raise exceptions.MotifStateTreeException(
                "cannot supply both a state and motif name to add a state to "
                "a motif state tree")

        if ms is None and m_name is None:
            raise exceptions.MotifStateTreeException(
                "must supply a motif state object or motif name to add_state")

        if ms is not None:
            for n in self.tree.nodes:
                if n.data.ref_state.uuid == ms.uuid:
                    raise exceptions.MotifStateTreeException(
                        "cannot add state: " + ms.name + " to tree as its uuid is " +
                        "already present in the tree")

    def _get_parent_node(self, parent_index):
        """
        gets node that serve as the parent of the current motif being added
        to the tree.

        :param parent_index: the tree index corresponding to the requested motif
        :type parent_index: int

        :return: tree node of parent
        :rtype: TreeNode
        """

        parent = self.tree.last_node

        if parent_index != -1:
            try:
                parent = self.tree.get_node(parent_index)
            except exceptions.TreeIndexException:
                raise exceptions.MotifStateTreeException(
                    "parent_index supplied: " + str(parent_index) + " does not " +
                    "exist in current motif state tree")

        return parent

    def _get_parent_available_ends(self, parent, parent_end_index,
                                   parent_end_name):
        """
        Gets the available ends of the parent that the current motif can align
        to. If either parent_end_index or parent_end_name are specified, checks
        to see if that position is available.

        :param parent: the TreeNode of parent
        :type parent: TreeNode

        :param parent_end_index: which end this motif will be aligned to on
            parent
        :type parent_end_index: int

        :param parent_end_name: the name instead of the index of the end the
            current motif will align to
        :type parent_end_name: str

        :return: list of available parent end indexes that meet the specified
            constraints
        :rtype: list of ints
        """

        if parent is None:
            return []

        if parent_end_index != -1 and parent_end_name is not None:
            raise exceptions.MotifStateTreeException(
                "cannot supply parent_end_index and parent_end_name together")

        elif parent_end_name is not None:
            try:
                parent_end_index = parent.data.get_end_index(name=parent_end_name)
            except:
                raise exceptions.MotifStateTreeException(
                    "cannot find parent_end_name: " + parent_end_name + " in "
                    "parent motif: " + parent.data.name())

            if parent_end_index == parent.data.block_end_add():
                raise exceptions.MotifStateTreeException(
                    "cannot add state: to tree as the parent_end_name" +
                    " supplied is blocked see class MotifState")

            available = parent.available_pos(parent_end_index)
            if not available:
                raise exceptions.MotifStateTreeException(
                    "cannot add state to tree as the end " +
                    "you are trying to add it to is already filled or does "
                    "not exist")

            return [parent_end_index]

        elif parent_end_index != -1:
            if parent_end_index == parent.data.block_end_add():
                raise exceptions.MotifStateTreeException(
                    "cannot add state: to tree as the parent_end_index" +
                    " supplied is blocked see class MotifState")

            available = parent.available_pos(parent_end_index)
            if not available:
                raise exceptions.MotifStateTreeException(
                    "cannot add state to tree as the end " +
                    "you are trying to add it to is already filled or does "
                    "not exist")
            return [parent_end_index]

        else:
            avail_pos = parent.available_children_pos()
            avail_pos.remove(0)
            return avail_pos

    def _get_state_from_manager(self, m_name, m_end_name):
        """
        helper function for add_motif should not be called directly. calls
        resource manager to get motif to be added to tree by the name of
        the motif.

        :param m_name: name of the motif to add to the tree
        :type m_name: str

        :param m_end_name: name of the basepair end of the motif to align by.
        :type m_end_name: str

        """

        try:
            if m_end_name is not None:
                state = rm.manager.get_state(name=m_name, end_name=m_end_name)
            else:
                state = rm.manager.get_state(name=m_name)
        except exceptions.ResourceManagerException as e:
            raise exceptions.MotifStateTreeException(
                "cannot add state to tree, state cannot be found in resource "
                "manager")

        return state

    def add_state(self, state=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None, m_name=None, m_end_name=None):

        self._validate_arguments_to_add_state(state, m_name)

        parent = self._get_parent_node(parent_index)

        if parent is None:
            n_data = NodeData(state)
            return self.tree.add_data(n_data, len(state.end_states), -1, -1)

        if state is None and m_name is not None:
            state = self._get_state_from_manager(m_name, m_end_name)

        avail_pos = self._get_parent_available_ends(parent, parent_end_index,
                                                    parent_end_name)

        for p in avail_pos:
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

        parent = self._get_parent_node(parent_index)
        parent_avail_ends = self._get_parent_available_ends(
                                parent, parent_end_index, parent_end_name)

        if len(parent_avail_ends) > 1 or len(parent_avail_ends) == 0:
            parent_end_index = -1
        else:
            parent_end_index = parent_avail_ends[0]

        index_dict = {}
        for i, n in enumerate(mst):
            if i == 0:
                j = self.add_state(n.data.ref_state,
                                   parent_index, parent_end_index)
            else:
                ind = index_dict[n.parent_index()]
                pei = n.parent_end_index()
                j = self.add_state(n.data.ref_state,
                                   parent_index=ind,
                                   parent_end_index=pei)
            if j == -1:
                raise exceptions.MotifStateTreeException(
                    "failed to add a state in add_mst to the current "
                    "motif tree it is likely a steric clash, consider "
                    "turning off sterics")

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

    def replace_state(self, i, new_state):
        n = self.get_node(i)
        if len(new_state.end_states) !=  len(n.data.ref_state.end_states):
            raise ValueError(
                "attempted to replace a state with a different number of ends")

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

    #REMOVE FUNCTIONS   #######################################################
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

    #TREE WRAPPER      ########################################################
    def get_node(self, i):
        return self.tree.get_node(i)

    def last_node(self):
        return self.tree.last_node

    def next_level(self):
        self.tree.level += 1

    def next(self):
        return self.tree.next()

    #MOTIF TREE WRAPPER      ##################################################
    def to_motif_tree(self):
        mt = motif_tree.MotifTree(sterics=self.option('sterics'))
        for i, n in enumerate(self.tree.nodes):
            m = rm.manager.get_motif(name=n.data.ref_state.name,
                                     end_name = n.data.ref_state.end_names[0])

            if i == 0:
                motif.align_motif(n.data.cur_state.end_states[0],
                                  m.ends[0],
                                  m)
                j = mt.add_motif(m)
            else:
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

    def write_pdbs(self, name="nodes"):
        self.to_motif_tree().write_pdbs(name)

    #MISC              ########################################################
    def _steric_clash(self, new_data):
        for n in self.tree.nodes[::-1]:
            for b1 in n.data.cur_state.beads:
                for b2 in new_data.cur_state.beads:
                    dist = util.distance(b1, b2)
                    if dist < self.clash_radius:
                        return 1
        return 0

    def topology_to_str(self):
        s = ""
        for n in self.tree.nodes:
            s += n.data.ref_state.name + "," + n.data.ref_state.end_ids[0] + "," + \
                 str(n.parent_index()) + "," + str(n.parent_end_index())  +  " "
        s += "|"
        for c in self.connections:
            s += c.to_str() + " "
        return s


class NodeData(object):
    def __init__(self, ref_state):
        self.ref_state = ref_state
        self.cur_state = ref_state.copy()

    def get_end_state(self, name=None, id=None):
        return self.cur_state.get_end_state(name, id)

    def get_end_index(self, name=None, id=None):
        return self.cur_state.get_end_index(name, id)

    def name(self):
        return self.cur_state.name

    def block_end_add(self):
        return self.cur_state.block_end_add


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














