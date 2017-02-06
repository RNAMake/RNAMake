import base
import option
import settings
import graph
import motif
import motif_graph
import resource_manager as rm
import exceptions
import util
import copy
from motif_state_node import NodeData


class MotifStateGraph(base.Base):
    def __init__(self, mg=None):
        super(self.__class__, self).__init__()
        self.setup_options_and_constraints()
        self.graph = graph.GraphStatic()
        self.clash_radius = settings.CLASH_RADIUS
        self.aligned = {}

        self.align_list = []
        self.update_align_list = 1

        if mg:
            self._setup_from_mg(mg)

    def setup_options_and_constraints(self):
        options = {'sterics': 1}

        self.options = option.Options(options)
        self.constraints = {}

    def _setup_from_mg(self, mg):
        self.aligned = copy.deepcopy(mg.aligned)
        max_index = 0
        for n in mg.graph.nodes:
            n_data = NodeData(n.data.get_state())
            self.graph.add_data(n_data, -1, -1, -1, len(n.data.ends),
                                orphan=1, index=n.index)

            if n.index > max_index:
                max_index = n.index

        self.graph.index = max_index+1

        for c in mg.graph.connections:
            self.graph.connect(c.node_1.index, c.node_2.index,
                               c.end_index_1, c.end_index_2)

        self.update_align_list = 1


    def __len__(self):
        return len(self.graph)

    def __iter__(self):
        self.graph.__iter__()
        return self

    def next(self):
        return self.graph.next()

    #ADD FUNCTIONS      #######################################################
    def _validate_arguments_to_add_state(self, ms, m_name):
        """
        makes sure the add_motif_state function is called correctly

        :param m: motif state to add to graph
        :type m: Motif object

        :param m_name: name of motif to add
        :type m_name: str

        :return: None
        """

        if ms is not None and m_name is not None:
            raise exceptions.MotifStateGraphException(
                "cannot supply both a state and motif name to add a state to "
                "a motif state graph")

        if ms is None and m_name is None:
            raise exceptions.MotifStateGraphException(
                "must supply a motif state object or motif name to add_state")

        if ms is not None:
            for n in self.graph.nodes:
                if n.data.ref_state.uuid == ms.uuid:
                    raise exceptions.MotifStateGraphException(
                        "cannot add state: " + ms.name + " to graph as its uuid is " +
                        "already present in the tree")

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
                ms = rm.manager.get_state(name=m_name, end_name=m_end_name)
            else:
                ms = rm.manager.get_state(name=m_name)
        except exceptions.ResourceManagerException as e:
            raise exceptions.MotifStateGraphException(
                "cannot add motif to graph, motif cannot be found in resource "
                "manager")

        return ms

    def _get_parent_node(self, parent_index):
        """
        gets node that serve as the parent of the current motif being added
        to the tree.

        :param parent_index: the tree index corresponding to the requested motif
        :type parent_index: int

        :return: tree node of parent
        :rtype: TreeNode
        """

        parent = self.graph.last_node

        if parent_index != -1:
            try:
                parent = self.graph.get_node(parent_index)
            except exceptions.GraphIndexException:
                raise exceptions.MotifStateGraphException(
                    "parent_index supplied: " + str(parent_index) + " does not " +
                    "exist in current motif graph")

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
            raise exceptions.MotifStateGraphException(
                "cannot supply parent_end_index and parent_end_name together")

        elif parent_end_name is not None:
            try:
                parent_end_index = parent.data.get_end_index(name=parent_end_name)
            except:
                raise exceptions.MotifStateGraphException(
                    "cannot find parent_end_name: " + parent_end_name + " in "
                    "parent motif: " + parent.data.name())

            if parent_end_index == parent.data.block_end_add():
                raise exceptions.MotifStateGraphException(
                    "cannot add state: to graph as the parent_end_name" +
                    " supplied is blocked see class MotifState")

            available = parent.available_pos(parent_end_index)
            if not available:
                raise exceptions.MotifStateGraphException(
                    "cannot add state to graph as the end " +
                    "you are trying to add it to is already filled or does "
                    "not exist")

            if self.connections.in_connection(parent.index, parent_end_name):
                raise exceptions.MotifStateGraphException(
                    "cannot add motif to graph as the end " +
                    "you are trying to add it to is in a connection")

            return [parent_end_index]

        elif parent_end_index != -1:
            if parent_end_index == parent.data.block_end_add():
                raise exceptions.MotifStateGraphException(
                    "cannot add state: to graph as the parent_end_index" +
                    " supplied is blocked see class MotifState")

            available = parent.available_pos(parent_end_index)
            if not available:
                raise exceptions.MotifStateGraphException(
                    "cannot add state to graph as the end " +
                    "you are trying to add it to is already filled or does "
                    "not exist")

            parent_end_name = parent.data.cur_state.end_names[parent_end_index]

            if self.connections.in_connection(parent.index, parent_end_name):
                raise exceptions.MotifStateGraphException(
                    "cannot add motif to graph as the end " +
                    "you are trying to add it to is in a connection")

            return [parent_end_index]

        else:
            avail_pos = parent.available_children_pos()

            final_avail_pos = []
            for p in avail_pos:
                if p == parent.data.block_end_add:
                    continue
                final_avail_pos.append(p)

            return final_avail_pos

    def add_state(self, ms=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None, m_name=None, m_end_name=None,
                  orphan=0):

        self._validate_arguments_to_add_state(ms, m_name)
        parent = self._get_parent_node(parent_index)

        if ms is None and m_name is not None:
            ms = self._get_state_from_manager(m_name, m_end_name)

        if parent is None or orphan:
            n_data = NodeData(ms)

            pos = self.graph.add_data(n_data, -1, -1, -1,
                                      len(ms.end_states), orphan=1)
            self.aligned[pos] = 0
            self.update_align_list = 1
            return pos


        avail_pos = self._get_parent_available_ends(parent, parent_end_index,
                                                    parent_end_name)


        for p in avail_pos:
            n_data = NodeData(ms)
            motif.get_aligned_motif_state(parent.data.cur_state.end_states[p],
                                          n_data.cur_state,
                                          n_data.ref_state)

            if self.option('sterics') and self._steric_clash(n_data):
                continue

            pos = self.graph.add_data(n_data, parent.index, p, 0, len(ms.end_states))
            self.aligned[pos] = 1
            self.update_align_list = 1
            return pos

        return -1

    def _get_connection_end(self, node, bp_name):
        node_end_index = -1

        if bp_name != "":
            ei = node.data.get_end_index(name=bp_name)
            if not node.available_pos(ei):
                raise exceptions.MotifStateGraphException(
                    "cannot add connection with " + str(node.index) + " and "
                    "end name " + bp_name + " as this end is not available")

            node_end_index = ei
        else:
            node_indexes = node.available_children_pos()

            if len(node_indexes) > 1:
                raise exceptions.MotifStateGraphException(
                    "cannot connect nodes " + str(node.index) + " its unclear "
                    " which ends to attach")

            if len(node_indexes) == 0:
                raise exceptions.MotifStateGraphException(
                    "cannot connect nodes " + str(node.index) + " there are "
                    "no ends free ends to attach too")

            node_end_index = node_indexes[0]
        return node_end_index

    def add_connection(self, i, j, i_bp_name="", j_bp_name=""):
        node_i = self.get_node(i)
        node_j = self.get_node(j)

        node_i_ei = self._get_connection_end(node_i, i_bp_name)
        node_j_ei = self._get_connection_end(node_j, j_bp_name)

        self.graph.connect(i, j, node_i_ei, node_j_ei)

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
            n.data.cur_state.uuid = uuid
            n.data.ref_state.uuid = uuid

        self._align_states(i)

    #REMOVE FUNCTIONS   #######################################################
    def remove_state(self, pos=-1):
        if pos == -1:
            pos = self.last_node().index
        n = self.graph.get_node(pos)
        self.graph.remove_node(pos)
        del self.aligned[pos]
        self.update_align_list = 1

    def remove_node_level(self, level=None):
        if level is None:
            level = self.graph.level

        r = range(1, len(self.graph.nodes))
        for i in r[::-1]:
            if self.graph.nodes[i].level >= level:
                self.remove_state(self.graph.nodes[i].index)

    #GRAPH WRAPPER      #######################################################
    def increase_level(self):
        self.graph.increase_level()

    def decrease_level(self):
        self.graph.decrease_level()

    def last_node(self):
        return self.graph.last_node

    def get_node(self, i=None, uuid=None, m_name=None):
        if i is not None:
            return self.graph.get_node(i)

        node = None
        for n in self.graph.nodes:
            if n.data.uuid() == uuid:
                return n
            if n.data.name() == m_name and m_name is not None:
                if node is not None:
                    raise exceptions.MotifStateGraphException(
                        "cannot get node with motif name: " + m_name + " there "
                        "are more then one")
                node = n

        if node is None:
            raise exceptions.MotifStateGraphException(
                "could not find node with uuid " + str(uuid) + " and m_name: "
                + str(m_name))

        return node

    #MOTIF GRAPH WRAPPER     ##################################################
    def to_motif_graph(self):
        align_list = self._get_align_list()
        non_aligned = self.get_not_aligned_nodes()
        mg = motif_graph.MotifGraph()
        mg.option('sterics', 0)
        seen_connections = {}
        index_hash = {}

        for n in align_list:
            m = rm.manager.get_motif(name=n.data.ref_state.name,
                                     end_name = n.data.ref_state.end_names[0])

            if n in non_aligned:
                motif.align_motif(n.data.cur_state.end_states[0],
                                  m.ends[0],
                                  m)
                j = mg.add_motif(m, orphan=1)
            else:
                c = n.connections[0]
                seen_connections[c] = 1
                parent = c.partner(n.index)
                parent_end_index = c.end_index(parent.index)
                j = mg.add_motif(m, parent_index=index_hash[parent.index],
                                 parent_end_index=parent_end_index)

            index_hash[n.index] = j


        for c in self.graph.connections:
            if c in seen_connections:
               continue
            mg.add_connection(index_hash[c.node_1.index],
                              index_hash[c.node_2.index],
                              c.node_1.data.end_name(c.end_index_1),
                              c.node_2.data.end_name(c.end_index_2))


        mg.update_indexes(index_hash)
        mg.option('sterics', self.option('sterics'))
        return mg

    def nodes_to_pdbs(self):
        self.to_motif_graph().nodes_to_pdbs()

    #MISC               #######################################################
    def _steric_clash(self, new_data):
        for n in self.graph.nodes[::-1]:
            for b1 in n.data.cur_state.beads:
                for b2 in new_data.cur_state.beads:
                    dist = util.distance(b1, b2)
                    if dist < self.clash_radius:
                        return 1
        return 0

    def _get_align_list(self):
        if not self.update_align_list:
            return self.align_list

        non_aligned_nodes = self.get_not_aligned_nodes()
        self.align_list = []
        used_nodes = {}
        for start in non_aligned_nodes:
            open = [start]
            seen_nodes = {}

            while open:
                n = open.pop(0)
                seen_nodes[n] = 1
                if n.index == start.index:
                    self.align_list.append(n)
                    used_nodes[n] = 1
                else:
                    if n.connections[0] is None:
                        continue
                    c = n.connections[0]
                    parent = c.partner(n.index)

                    if parent not in used_nodes:
                        continue

                    self.align_list.append(n)

                used_nodes[n] = 1
                for i, c in enumerate(n.connections):
                    if i == n.data.block_end_add() or c is None:
                        continue

                    partner_n = c.partner(n.index)
                    if partner_n in seen_nodes or partner_n in used_nodes:
                        continue

                    # if something goes wrong check this!
                    if c.end_index(partner_n.index) == partner_n.data.block_end_add():
                        open.append(partner_n)
                    elif len(partner_n.data.cur_state.end_states) == 1:
                        open.append(partner_n)

        self.update_align_list = 0
        return self.align_list

    def get_not_aligned_nodes(self):
        not_aligned = []
        for n in self.graph:
            if self.aligned[n.index] == 0:
                not_aligned.append(n)
        return not_aligned

    def _align_states(self, pos=-1):
        non_aligned_nodes = self.get_not_aligned_nodes()
        align_list = self._get_align_list()

        start = 1
        if pos != -1:
            start = 0

        for n in align_list:
            if start == 0:
                if n.index == pos:
                    start = 1
                else:
                    continue

            if n in non_aligned_nodes:
                continue

            parent = n.connections[0].partner(n.index)
            pei = n.connections[0].end_index(parent.index)

            motif.get_aligned_motif_state(parent.data.cur_state.end_states[pei],
                                          n.data.cur_state,
                                          n.data.ref_state)

