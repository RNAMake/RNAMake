from collections import namedtuple, defaultdict

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
from motif_state_node import NodeData


class MotifStateTree(base.Base):

    class _MotifStateTreePrinter(object):
        """
        A small private class to handle pretty printing a tree for visual
        inspection of the connectivity of the tree

        :param mst: the MotifStateTree instance to print out
        :type mst: MotifStateTree

        :attributes:

        `mst` : MotifStateTree
            The MotifStateTree instance to print out
        `levels` : Dict of key and values of ints
            Keeps track of which nodes are at each level in the tree. i.e.
            how many levels of seperation between the first node and the current
            node
        `node_pos` : Dict of key and values of ints
            Keeps track of the horizontal position of each node on the screen.
            A value of 10 would be 10 spaces from the left.
        `branch_length` : int
            The initial spacing between nodes from the same parent. For example
            if a parent had a node_pos of 100 and the branch_length was 25.
            Then child one would be at pos 75 and the other would be at 125.
        `start_pos` : int
            The horizontal position of the first node on the screen.
        `node_per_level` : Dict of key and values of ints
            Keeps track of how many nodes inhabit each level
        """

        def __init__(self, mst):
            self.mst = mst
            self.levels = {}
            self.node_pos = {}
            self.branch_length = 25
            self.start_pos = 100
            self.nodes_per_level = {}

            self._setup_node_positions()

        def _assign_node_levels(self):
            """
            calculates and stores the node level of each node in the mt. The
            head node has a level of 1 and its children have a level of 2 and
            so on.
            """
            for n in self.mst:
                if len(self.levels) == 0:
                    self.levels[n.index] = 1
                else:
                    parent_level = self.levels[n.parent_index()]
                    self.levels[n.index] = parent_level + 1

        def _setup_node_positions(self):
            """
            Setups up the horizontal position of each node on the screen.
            Position is prograted from the start_pos with the first node. If
            the node only has one child, that child retains the same position,
            but if it has two the children will be seperated by 2* the branch
            length.
            """

            self._assign_node_levels()

            self.nodes_per_level = defaultdict(int)
            for n in self.mst:
                self.nodes_per_level[self.levels[n.index]] += 1

            for i, n in enumerate(self.mst):
                if i == 0:
                    self.node_pos[n.index] = self.start_pos

                children = []
                for c in n.children:
                    if c is not None:
                        children.append(c)
                if len(children) == 1:
                    self.node_pos[children[0].index] = self.node_pos[n.index]
                elif len(children) == 2:
                    level = self.levels[n.index]
                    nodes_per_level = self.nodes_per_level[level+1]
                    extra = nodes_per_level - 2
                    parent_pos = self.node_pos[n.index]
                    if extra == 0:
                        self.node_pos[children[0].index] = parent_pos - self.branch_length
                        self.node_pos[children[1].index] = parent_pos + self.branch_length
                    else:
                        self.node_pos[children[0].index] = parent_pos - self.branch_length / extra
                        self.node_pos[children[1].index] = parent_pos + self.branch_length / extra
                elif len(children) > 2:
                    raise exceptions.MotifTreeException(
                        "Greater then two children is not supported for pretty_printing")

        def _print_level(self, nodes):
            """
            creates a string of formatted information of each node in on the
            current level.

            :param nodes: the nodes on the current level
            :type nodes: list of Tree nodes
            """
            nodes_and_pos = []
            for n in nodes:
                nodes_and_pos.append([n, self.node_pos[n.index]])

            s = ""
            nodes_and_pos.sort(key=lambda x: x[1])

            strings  = []
            for n, pos in nodes_and_pos:
                strs = []
                if n.parent is not None:
                    parent_end_index = n.parent_end_index()
                    parent_end_name = n.parent.data.end_name(parent_end_index)
                    strs.append("|")
                    strs.append("E" + str(parent_end_index) +  " - " + \
                                parent_end_name)
                    strs.append("|")

                strs.append("N" + str(n.index) + " - " + n.data.name())
                strs.append("|  - " + n.data.end_name(0))
                strings.append(strs)

            transposed_strings = []
            for i in range(len(strings[0])):
                transposed = []
                for strs in strings:
                    transposed.append(strs[i])
                transposed_strings.append(transposed)


            for strs in transposed_strings:
                current_pos = 0
                j = 0
                for n, pos in nodes_and_pos:
                    diff = pos - current_pos
                    cur_s = '%'+str(diff)+'s'
                    s += cur_s % ("")
                    s += strs[j]
                    current_pos = pos + len(strs[j])
                    j += 1

                s+= "\n"

            current_pos = 0
            hit = 0
            for n, pos in nodes_and_pos:
                children = []
                for c in n.children:
                    if c is not None:
                        children.append(c)
                if len(children) > 1:
                    hit = 1
                    min = self.node_pos[children[0].index]
                    max = self.node_pos[children[-1].index]

                    diff = min+1 - current_pos
                    cur_s = '%'+str(diff)+'s'
                    s += cur_s % ("")
                    for i in range(max-min-1):
                        s += "_"
                    current_pos = max

            if hit:
                s += "\n"

            return s

        def print_tree(self):
            """
            Actually generates the formatted string for the entire tree. This
            is the only function that should be called in the normal use of
            this class

            :returns: formatted string of entire tree
            :rtype: str
            """

            nodes_per_level = defaultdict(list)
            for n in self.mst:
                nodes_per_level[self.levels[n.index]].append(n)

            found = 1
            level = 1
            s = "\n"
            while found:
                if level not in nodes_per_level:
                    break

                node_level = nodes_per_level[level]
                s += self._print_level(node_level)

                level += 1
            return s


    #SETUP FUNCTIONS ##########################################################
    def __init__(self, mt=None, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.tree = tree.TreeStatic()
        self.clash_radius = settings.CLASH_RADIUS
        self.connections = motif_connection.MotifConnections()

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

        self.connections = mt.connections.copy()

    def copy(self):
        mst = MotifStateTree()
        new_tree = self.tree.copy()
        mst.tree = new_tree
        mst.connections = self.connections.copy()
        return mst
        #return str_to_motif_state_tree(self.topology_to_str(), sterics=0)

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

            if self.connections.in_connection(parent.index, parent_end_name):
                raise exceptions.MotifStateTreeException(
                    "cannot add motif to tree as the end " +
                    "you are trying to add it to is in a connection")

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

            parent_end_name = parent.data.cur_state.end_names[parent_end_index]

            if self.connections.in_connection(parent.index, parent_end_name):
                raise exceptions.MotifStateTreeException(
                    "cannot add motif to tree as the end " +
                    "you are trying to add it to is in a connection")

            return [parent_end_index]

        else:
            avail_pos = parent.available_children_pos()
            avail_pos.remove(0)

            final_avail_pos = []
            for p in avail_pos:
                pen = parent.data.cur_state.end_names[p]
                if self.connections.in_connection(parent.index, pen):
                    continue
                if p == parent.data.block_end_add:
                    continue
                final_avail_pos.append(p)

            return final_avail_pos

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

    def _get_connection_end(self, node, bp_name):

        node_end_index = -1

        if bp_name != "":
            ei = node.data.get_end_index(name=bp_name)
            if not node.available_pos(ei):
                raise exceptions.MotifStateTreeException(
                    "cannot add connection with " + str(node.index) + " and "
                    "end name " + bp_name + " as this end is not available")

            if self.connections.in_connection(node.index, bp_name):
                raise exceptions.MotifStateTreeException(
                    "cannot add connection with " + node.index + " and end "
                    "name " + bp_name + " as this end is already in a "
                    "connection")

            node_end_index = ei
        else:
            node_indexes = node.available_children_pos()
            node_indexes.remove(0)

            if len(node_indexes) > 1:
                raise exceptions.MotifStateTreeException(
                    "cannot connect nodes " + str(node.index) + " its unclear "
                    " which ends to attach")

            if len(node_indexes) == 0:
                raise exceptions.MotifStateTreeException(
                    "cannot connect nodes " + str(node.index) + " there are "
                    "no ends free ends to attach too")

            node_index_name = node.data.cur_state.end_names[node_indexes[0]]
            if self.connections.in_connection(node.index, node_index_name):
                raise exceptions.MotifStateTreeException(
                    "cannot add connection with " + str(node.index) + " and end "
                    "name " + node_index_name + " as this end is already in a "
                    "connection")

            node_end_index = node_indexes[0]

        return node_end_index

    def add_connection(self, i, j, i_bp_name="", j_bp_name=""):
        node_i = self.get_node(i)
        node_j = self.get_node(j)

        node_i_ei = self._get_connection_end(node_i, i_bp_name)
        node_j_ei = self._get_connection_end(node_j, j_bp_name)

        node_i_end_name = node_i.data.cur_state.end_names[node_i_ei]
        node_j_end_name = node_j.data.cur_state.end_names[node_j_ei]

        self.connections.add_connection(i, j, node_i_end_name, node_j_end_name)

    def replace_state(self, i, new_state, keep_uuid=1):
        n = self.get_node(i)
        if len(new_state.end_states) !=  len(n.data.ref_state.end_states):
            raise ValueError(
                "attempted to replace a state with a different number of ends")

        old_state = n.data.ref_state
        if keep_uuid:
            uuid = n.data.ref_state.uuid

        n.data.ref_state = new_state
        n.data.cur_state = new_state.copy()

        if keep_uuid:
            n.data.ref_state.uuid = uuid
            n.data.cur_state.uuid = uuid

        for i in range(len(old_state.end_names)):
            if self.connections.in_connection(n.index, old_state.end_names[i]):
                self.connections.update_connection_name(n.index,
                                                        old_state.end_names[i],
                                                        new_state.end_names[i])


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
    def get_node(self, i=None, uuid=None):
        if i is not None:
            return self.tree.get_node(i)
        elif uuid is not None:
            for n in self.tree:
                if n.data.cur_state.uuid == uuid:
                    return n
            raise exceptions.MotifStateTreeException(
                "cannot find motif state with uuid")

    def last_node(self):
        return self.tree.last_node

    def increase_level(self):
        self.tree.increase_level()

    def decrease_level(self):
        self.tree.decrease_level()

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

    def get_structure(self):
        mt = self.to_motif_tree()
        return mt.get_structure()

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

    def to_pretty_str(self):
        printer = self._MotifStateTreePrinter(self)
        return printer.print_tree()




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














