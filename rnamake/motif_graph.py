import base
import option
import settings
import motif
import util
import residue
import motif_type
import graph
import resource_manager as rm
import motif_merger
import copy
import exceptions
import steric_lookup

from collections import defaultdict

class MotifGraph(base.Base):
    class _MotifGraphBuildPoint(object):
        def __init__(self, node, end_index):
            self.node = node
            self.end_index = end_index


    class _MotifGraphPrinter(object):
        """
        A small private class to handle pretty printing a tree for visual
        inspection of the connectivity of the tree

        :param mt: the MotifTree instance to print out
        :type mt: MotifTree

        :attributes:

        `mt` : MotifTree
            The MotifTree instance to print out
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

        def __init__(self, mg):
            self.mg = mg
            self.nodes = []
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

            nodes = self.mg.get_not_aligned_nodes()

            if len(nodes) == 0:
                raise exceptions.MotifGraphException(
                  "cannot find a place to start printing in motif_graph"
                  " to_pretty_str")

            start = nodes[0].index
            i = 0
            for n in graph.transverse_graph(self.mg.graph, start, directed=0):
                if len(self.levels) == 0:
                    self.levels[n.index] = 1
                else:
                    c = n.connections[0]
                    parent = c.partner(n.index)
                    self.levels[n.index] = self.levels[parent.index] + 1
                self.nodes.append(n)

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
            for n in self.nodes:
                self.nodes_per_level[self.levels[n.index]] += 1

            for i, n in enumerate(self.nodes):
                if i == 0:
                    self.node_pos[n.index] = self.start_pos

                children = []
                for j, c in enumerate(n.connections):
                    if j == 0:
                        continue
                    if c is not None:
                        child = c.partner(n.index)
                        if child.connections[0] != c:
                            continue

                        if child.index in self.node_pos:
                            continue
                        else:
                            children.append(child)

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

        def _print_pos(self):
            """
            used for testing purposes, prints out the position of where each
            node will appear
            """

            nodes_per_level = defaultdict(list)
            for n in self.mt:
                nodes_per_level[self.levels[n.index]].append(n)

            level = 1
            found = 1
            s = "\n"
            while found:
                if level not in nodes_per_level:
                    break

                node_level = nodes_per_level[level]
                nodes_and_pos = []
                for n in node_level:
                    nodes_and_pos.append([n, self.node_pos[n.index]])

                nodes_and_pos.sort(key=lambda x: x[1])
                current_pos = 0
                for n, pos in nodes_and_pos:
                    diff = pos - current_pos
                    cur_s = '%'+str(diff)+'s'
                    s += cur_s % (n.index)
                    current_pos = pos
                s += "\n"

                level += 1
            return s

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
                if n.parent() is not None:
                    parent_end_index = n.parent_end_index()
                    parent_end_name = n.parent().data.ends[parent_end_index].name()
                    strs.append("|")
                    strs.append("E" + str(parent_end_index) +  " - " + \
                                parent_end_name)
                    strs.append("|")

                strs.append("N" + str(n.index) + " - " + n.data.name)
                strs.append("|  - " + n.data.ends[0].name())
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
                for i, c in enumerate(n.connections):
                    if i == 0:
                        continue
                    if c is not None:
                        children.append(c.partner(n.index))
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

        def print_graph(self):
            """
            Actually generates the formatted string for the entire tree. This
            is the only function that should be called in the normal use of
            this class

            :returns: formatted string of entire tree
            :rtype: str
            """

            nodes_per_level = defaultdict(list)
            for n in self.nodes:
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
    def __init__(self, mg_str="", top_str=""):
        super(self.__class__, self).__init__()
        self.setup_options_and_constraints()
        self.graph = graph.GraphStatic()
        self.clash_radius = settings.CLASH_RADIUS
        self.merger = motif_merger.MotifMerger()
        self.aligned = {}

        if top_str != "":
            self._setup_from_top_str(top_str)

        if mg_str != "":
            self._setup_from_str(mg_str)

    def __len__(self):
        return len(self.graph)

    def __iter__(self):
        self.graph.__iter__()
        return self

    def next(self):
        return self.graph.next()

    def _setup_from_top_str(self, s):
        self.option('sterics', 0)
        spl = s.split("&")
        node_spl = spl[0].split("|")
        max_index = 0
        for i, n_spl in enumerate(node_spl[:-1]):
            sspl = n_spl.split(",")
            if rm.manager.contains_motif(name=sspl[0], end_name=sspl[1]):
                m = rm.manager.get_motif(name=sspl[0], end_name=sspl[1])
            else:
                raise exceptions.MotifGraphException(
                    "Unknown motif name: " + sspl[0] + " end_name: " + sspl[1] +\
                    " did you forget to add your custom motifs to the resource "
                    "manager ")

            m_copy = m.copy()
            m_copy.get_beads(m_copy.ends)
            pos = self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends),
                                      orphan=1, index=int(sspl[2]))
            self.aligned[int(sspl[2])] = int(sspl[3])

            if int(sspl[2]) > max_index:
                max_index = int(sspl[2])

        self.graph.index = max_index+1

        con_spl = spl[1].split("|")
        for c_str in con_spl[:-1]:
            c_spl = c_str.split(",")
            self.graph.connect(int(c_spl[0]), int(c_spl[1]),
                               int(c_spl[2]), int(c_spl[3]))

        start = -1
        for k,v in self.aligned.iteritems():
            if v == 0:
                start = k

        if start == -1:
            raise ValueError("cannot find a place to start in rebuilding motif_graph"
                             " from string")


        for n in graph.transverse_graph(self.graph, start):
            if n.index == start:
                self.merger.add_motif(n.data)
                continue
            if n.connections[0] is None:
                continue

            c = n.connections[0]
            parent = c.partner(n.index)
            parent_end_index = c.end_index(parent.index)

            m_added = motif.get_aligned_motif(parent.data.ends[parent_end_index],
                                              n.data.ends[0],
                                              n.data)
            n.data = m_added

            self.merger.add_motif(n.data, n.data.ends[0],
                                  parent.data, parent.data.ends[parent_end_index])

    def _setup_from_str(self, s):
        self.option('sterics', 0)
        spl = s.split("FAF")
        node_spl = spl[0].split("KAK")
        max_index = 0
        for i, n_str in enumerate(node_spl[:-1]):
            n_spl = n_str.split("^")
            m = motif.str_to_motif(n_spl[0])
            if len(m.ends) > 0:
                m.get_beads([m.ends[0]])
            else:
                m.get_beads()
            self.graph.add_data(m, -1, -1, -1, len(m.ends),
                                orphan=1, index=int(n_spl[1]))
            self.aligned[int(n_spl[1])] = int(n_spl[2])

            if int(n_spl[1]) > max_index:
                max_index = int(n_spl[1])

        self.graph.index = max_index+1

        con_spl = spl[1].split("|")
        for c_str in con_spl[:-1]:
            c_spl = c_str.split(",")
            self.graph.connect(int(c_spl[0]), int(c_spl[1]),
                               int(c_spl[2]), int(c_spl[3]))

        start = -1
        for k,v in self.aligned.iteritems():
            if v == 0:
                start = k

        if start == -1:
            raise ValueError("cannot find a place to start in rebuilding motif_graph"
                             " from string")

        for n in graph.transverse_graph(self.graph, start):
            if n.index == start:
                self.merger.add_motif(n.data)
                continue
            if n.connections[0] is None:
                continue

            c = n.connections[0]
            parent = c.partner(n.index)
            parent_end_index = c.end_index(parent.index)

            self.merger.add_motif(n.data, n.data.ends[0],
                                  parent.data, parent.data.ends[parent_end_index])

    def setup_options_and_constraints(self):
        options = {'sterics': 1}

        self.options = option.Options(options)
        self.constraints = {}

    def copy(self):
        mg = MotifGraph()
        new_graph = self.graph.copy()
        mg.graph = new_graph
        mg.merger = self.merger.copy([n.data for n in new_graph.nodes])
        mg.aligned = copy.deepcopy(self.aligned)
        return mg

    #ADD FUNCTIONS      #######################################################
    def _validate_arguments_to_add_motif(self, m, m_name):
        """
        makes sure the add_motif function is called correctly

        :param m: motif to add to graph
        :type m: Motif object

        :param m_name: name of motif to add
        :type m_name: str

        :return: None
        """

        if m is not None and m_name is not None:
            raise exceptions.MotifGraphException(
                "cannot supply both a motif and motif name to add a motif to "
                "a motif graph")

        if m is None and m_name is None:
            raise exceptions.MotifGraphException(
                "must supply a motif object or motif name to add_motif")

        if m is not None:
            for n in self.graph.nodes:
                if n.data.id == m.id:
                    raise exceptions.MotifGraphException(
                        "cannot add motif: " + m.name + " to graph as its uuid is " +
                        "already present in the graph")

    def _get_motif_from_manager(self, m_name, m_end_name):
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
                m = rm.manager.get_motif(name=m_name, end_name=m_end_name)
            else:
                m = rm.manager.get_motif(name=m_name)
        except exceptions.ResourceManagerException as e:
            raise exceptions.MotifGraphException(
                "cannot add motif to graph, motif cannot be found in resource "
                "manager")

        return m

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
                raise exceptions.MotifGraphException(
                    "parent_index supplied: " + str(parent_index) + " does not " +
                    "exist in current motif graph")

        return parent

    def _get_parent_available_ends(self, parent, parent_end_index,
                                   parent_end_name):
        """
        Gets the available ends of the parent that the current motif can align
        to. If either parent_end_index or parent_end_name are specified, checks
        to see if that position is available.

        :param parent: the GraphNode of parent
        :type parent: GraphNode

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
            raise exceptions.MotifGraphException(
                "cannot supply parent_end_index and parent_end_name together")

        elif parent_end_name is not None:
            parent_ends = parent.data.get_basepair(name=parent_end_name)
            if len(parent_ends) == 0:
                raise exceptions.MotifGraphException(
                    "cannot find parent_end_name: " + parent_end_name + " in "
                    "parent motif: " + parent.data.name)
            if len(parent_ends) > 1:
                raise exceptions.MotifGraphException(
                    "more then one end was found with parent_end_name: " +
                    parent_end_name + " in parent motif: " + parent.data.name)

            parent_end = parent_ends[0]
            parent_end_index = parent.data.ends.index(parent_end)

            if parent_end_index == parent.data.block_end_add:
                raise exceptions.MotifGraphException(
                    "cannot add motif: to graph as the parent_end_name" +
                    " supplied is blocked see class Motif")

            available = parent.available_pos(parent_end_index)
            if not available:
                raise exceptions.MotifGraphException(
                    "cannot add motif to tree as the end " +
                    "you are trying to add it to is already filled or does "
                    "not exist")

            return [parent_end_index]

        elif parent_end_index != -1:
            if parent_end_index == parent.data.block_end_add:
                raise exceptions.MotifGraphException(
                    "cannot add motif: to tree as the parent_end_index" +
                    " supplied is blocked see class Motif")

            available = parent.available_pos(parent_end_index)
            if not available:
                raise exceptions.MotifGraphException(
                    "cannot add motif to tree as the end " +
                    "you are trying to add it to is already filled or does "
                    "not exist")

            return [parent_end_index]

        else:
            avail_pos = parent.available_children_pos()

            final_avail_pos = []
            for p in avail_pos:
                if p == parent.data.block_end_add:
                    continue
                final_avail_pos.append(p)

            return final_avail_pos

    def _add_motif_to_graph(self, m, parent, parent_end_index):
        m.new_res_uuids()

        if parent is None:
            pos = self.graph.add_data(m, -1, -1, -1, len(m.ends), orphan=1)
            if pos != -1:
                self.merger.add_motif(m)
                self.aligned[pos] = 0
            return pos

        else:
            pos = self.graph.add_data(m, parent.index, parent_end_index,
                                      0, len(m.ends))

            if pos != -1:
                self.merger.add_motif(m, m.ends[0],
                                      parent.data, parent.data.ends[parent_end_index])
                self.aligned[pos] = 1
            return pos

    def add_motif(self, m=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None, m_name=None, m_end_name=None,
                  orphan=0):

        self._validate_arguments_to_add_motif(m, m_name)
        parent = self._get_parent_node(parent_index)

        if m is None and m_name is not None:
            m = self._get_motif_from_manager(m_name, m_end_name)
        else:
            if not rm.manager.contains_motif(name=m.name,
                                             end_name=m.ends[0].name()):
                rm.manager.register_motif(m)

        if parent is None or orphan:
            m_copy = m.copy()
            m_copy.get_beads([m_copy.ends[0]])

            return self._add_motif_to_graph(m_copy, None, None)

        avail_pos = self._get_parent_available_ends(parent, parent_end_index,
                                                    parent_end_name)

        for p in avail_pos:
            m_added = motif.get_aligned_motif(parent.data.ends[p], m.ends[0], m)
            if self.option('sterics'):
                if self._steric_clash(m_added):
                    continue

            return self._add_motif_to_graph(m_added, parent, p)

        return -1

    def add_motif_tree(self, mt, parent_index=-1, parent_end_name=""):

        if parent_index != -1:
            parent = self.get_node(parent_index)
            bps = parent.data.get_basepair(name=parent_end_name)
            if len(bps) == 0:
                raise ValueError("cannot find parent end in add_motif_tree")
            pei = parent.data.ends.index(bps[0])
        else:
            pei = -1

        index_hash = {}
        for i, n in enumerate(mt):
            m = rm.manager.get_motif(name=n.data.name, end_name=n.data.ends[0].name())
            if i == 0:
                j = self.add_motif(m, parent_index, pei)
            else:
                pi = index_hash[n.parent_index()]
                j = self.add_motif(m, pi, n.parent_end_index())
                if j == -1:
                    self.write_pdbs()
                    print m.name, self.option('sterics')
                    raise ValueError("cannot add_motif_tree")

            index_hash[n.index] = j

    def _get_connection_end(self, node, bp_name):
        node_end_index = -1

        if bp_name != "":
            ei = node.data.get_end_index(bp_name)
            if ei == node.data.block_end_add:
                raise exceptions.MotifGraphException(
                    "cannot add connection with " + str(node.index) + " and "
                    "end name " + bp_name + " as the end is blocked")

            if not node.available_pos(ei):
                raise exceptions.MotifGraphException(
                    "cannot add connection with " + str(node.index) + " and "
                    "end name " + bp_name + " as this end is not available")

            node_end_index = ei
        else:
            node_indexes = node.available_children_pos()

            if len(node_indexes) > 1:
                raise exceptions.MotifGraphException(
                    "cannot connect nodes " + str(node.index) + " its unclear "
                    " which ends to attach")

            if len(node_indexes) == 0:
                raise exceptions.MotifGraphException(
                    "cannot connect nodes " + str(node.index) + " there are "
                    "no ends free ends to attach too")

            node_end_index = node_indexes[0]
        return node_end_index

    def add_connection(self, i, j, i_bp_name="", j_bp_name=""):
        node_i = self.get_node(i)
        node_j = self.get_node(j)

        node_i_ei = self._get_connection_end(node_i, i_bp_name)
        node_j_ei = self._get_connection_end(node_j, j_bp_name)

        node_i_end_name = node_i.data.ends[node_i_ei].name()
        node_j_end_name = node_j.data.ends[node_j_ei].name()

        self.graph.connect(i, j, node_i_ei, node_j_ei)
        self.merger.connect_motifs(node_i.data, node_j.data,
                                   node_i.data.ends[node_i_ei],
                                   node_j.data.ends[node_j_ei])

    #REMOVE FUNCTIONS   #######################################################
    def remove_motif(self, pos):
        n = self.graph.get_node(pos)
        self.merger.remove_motif(n.data)
        self.graph.remove_node(pos)
        del self.aligned[pos]

    def remove_node_level(self, level=None):
        if level is None:
            level = self.graph.level

        r = range(1, len(self.graph.nodes))
        for i in r[::-1]:
            if self.graph.nodes[i].level >= level:
                self.remove_motif(self.graph.nodes[i].index)

    #DESIGNING          #######################################################
    def designable_secondary_structure(self):
        ss = self.merger.secondary_structure()

        for n in self.graph.nodes:
            if n.data.name != "HELIX.IDEAL":
                continue
            for r in n.data.residues():
                r_ss = ss.get_residue(uuid=r.uuid)
                if r_ss is not None:
                    r_ss.name = "N"

        return ss

    def replace_ideal_helices(self):
        found = 1
        while found:
            found = 0
            for n in self.graph.nodes:
                if n.data.mtype != motif_type.HELIX:
                    continue
                if len(n.data.residues()) == 4:
                    continue

                found = 1

                parent = None
                parent_end_index = None
                other = None
                other_end_index = None
                if n.connections[0] is not None:
                    parent = n.connections[0].partner(n.index)
                    parent_end_index = n.connections[0].end_index(parent.index)
                if n.connections[1] is not None:
                    other = n.connections[1].partner(n.index)
                    other_end_index = n.connections[1].end_index(other.index)

                name_spl = n.data.name.split(".")
                if len(name_spl) == 3:
                    count = int(name_spl[2])
                else:
                    count = 1
                i = n.index
                self.remove_motif(i)

                h = rm.manager.get_motif(name="HELIX.IDEAL")
                if parent is None:
                    h.get_beads([h.ends[0]])
                    pos = self._add_motif_to_graph(h, None, None)
                else:
                    m_added = motif.get_aligned_motif(parent.data.ends[parent_end_index],
                                                      h.ends[0], h)
                    pos = self._add_motif_to_graph(m_added, parent, parent_end_index)


                for j in range(0, count):
                    h = rm.manager.get_motif(name="HELIX.IDEAL")
                    parent =  self.graph.get_node(pos)
                    m_added = motif.get_aligned_motif(parent.data.ends[1], h.ends[0], h)
                    pos = self._add_motif_to_graph(m_added, parent, 1)

                if other:
                    self.graph.connect(pos, other.index, 1, other_end_index)
                    node = self.graph.get_node(pos)
                    self.merger.connect_motifs(node.data, other.data,
                                              node.data.ends[1],
                                             other.data.ends[other_end_index])
                    self.aligned[other.index] = 1

                break

    def replace_helix_sequence(self, ss):

        for n in self.graph.nodes:
            if n.data.mtype != motif_type.HELIX:
                continue
            ss_m = ss.motif(n.data.id)
            spl = ss_m.end_ids[0].split("_")
            new_name = spl[0][0] + spl[2][1] + "=" + spl[0][1] + spl[2][0]
            if new_name == n.data.name:
                continue

            m = rm.manager.get_motif(name=new_name)
            m.id = n.data.id
            org_res = n.data.residues()
            new_res = m.residues()

            for i in range(len(org_res)):
                new_res[i].uuid = org_res[i].uuid

            for i in range(len(n.data.basepairs)):
                m.basepairs[i].uuid = n.data.basepairs[i].uuid

            n.data = m

        self._align_motifs_all_motifs()

    def replace_motifs(self, motifs):
        for i, m in motifs.iteritems():
            self.merger.replace_motif(self.get_node(i).data, m)
            self.get_node(i).data = m

        self._align_motifs_all_motifs()

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
            if n.data.id == uuid:
                return n
            if n.data.name == m_name:
                if node is not None:
                    raise exceptions.MotifGraphException(
                        "cannot get node with motif name: " + m_name + " there "
                        "are more then one")
                node = n
        if node is None:
            raise exceptions.MotifGraphException(
                "could not find node with uuid " + str(uuid) + " and m_name: "
                + m_name)

        return node

    #MERGER WRAPPER     #######################################################
    def secondary_structure(self):
        return self.merger.secondary_structure()

    def get_structure(self):
        return self.merger.get_structure()

    #OUTPUTING          #######################################################
    def write_pdbs(self, name="node"):
        for n in self.graph.nodes:
            n.data.to_pdb(name + "." + str(n.index) + ".pdb")

    def to_pdb(self, name="test.pdb", renumber=-1, close_chain=0):
        return self.merger.get_structure().to_pdb(name, renumber=renumber,
                                                  close_chain=close_chain)

    def topology_to_str(self):
        s = ""
        con_str = ""
        seen_connections = {}
        for n in self.graph.nodes:
            s += n.data.name + "," + n.data.ends[0].name() + "," + str(n.index) + ","
            s += str(self.aligned[n.index]) + "|"
            for c in n.connections:
                if c is None:
                    continue
                key1 = str(c.node_1.index) + " " + str(c.node_2.index)
                key2 = str(c.node_2.index) + " " + str(c.node_1.index)
                if key1 in seen_connections or key2 in seen_connections:
                    continue
                seen_connections[key1] = 1
                con_str += str(c.node_1.index) + "," + str(c.node_2.index) + ","
                con_str += str(c.end_index_1) + "," + str(c.end_index_2) + "|"
        s += "&"
        s += con_str
        return s

    def to_str(self):
        s = ""
        con_str = ""
        seen_connections = {}
        for n in self.graph.nodes:
            s += n.data.to_str() + "^" + str(n.index) + "^"
            s += str(self.aligned[n.index]) + " KAK "
            for c in n.connections:
                if c is None:
                    continue
                key1 = str(c.node_1.index) + " " + str(c.node_2.index)
                key2 = str(c.node_2.index) + " " + str(c.node_1.index)
                if key1 in seen_connections or key2 in seen_connections:
                    continue
                seen_connections[key1] = 1
                con_str += str(c.node_1.index) + "," + str(c.node_2.index) + ","
                con_str += str(c.end_index_1) + "," + str(c.end_index_2) + "|"
        s += " FAF "
        s += con_str
        return s

    def to_pretty_str(self):
        printer = self._MotifGraphPrinter(self)
        return printer.print_graph()

    #GETTERS            #######################################################
    def get_build_points(self):
        build_points = []
        for n in self.graph.nodes:
            for i, c in enumerate(n.connections):
                if c is None and i != 0:
                    build_points.append(self._MotifGraphBuildPoint(n, i))

        return build_points

    def leafs_and_ends(self):
        build_points = []
        for n in self.graph.nodes:
             for i, c in enumerate(n.connections):
                if c is None and i != 0:
                    build_points.append([n, i])
        return build_points

    def get_beads(self, exclude_phos=1):
        pass

    def get_steric_lookup_table(self, exclude_phos=1):
        pass

    def get_end(self, pos=-1, m_name="", m_end_name=""):
        n = None
        if pos != -1:
            n = self.graph.get_node(pos)
        elif m_name != "":
            nodes = []
            for n in self.graph.nodes:
                if n.data.name == m_name:
                    nodes.append(n)
            if len(nodes) > 1:
                raise ValueError("cannot get end, too many motifs match name given "
                                 + m_name)
            n = nodes[0]

        if m_end_name != "":
            bps = n.data.get_basepair(name=m_end_name)
            if len(bps) == 0:
                raise ValueError("found motif but " + m_end_name + "is not is a "
                                                                   "basepair in it")
            end_index = n.data.get_end_index(name=m_end_name)
            end = n.data.ends[end_index]
            # check to see if the position is available
            self.graph.get_availiable_pos(n, end_index)
            return end

        else:
            avail_pos = n.available_children_pos()
            if len(avail_pos) > 1:
                raise ValueError("too many free ends to pick one in get_end")
            if len(avail_pos) == 0:
                raise ValueError("no ends available in get_end")

            end = n.data.ends[avail_pos[0]]
            return end

    def get_node_num(self, m_name):
        for n in self.graph.nodes:
            if n.data.name == m_name:
                return n.index
        return -1

    def get_not_aligned_nodes(self):
        not_aligned = []
        for n in self.graph:
            if self.aligned[n.index] == 0:
                not_aligned.append(n)
        return not_aligned

    def get_node_by_id(self, uuid):
        for n in self.graph.nodes:
            if n.data.id == uuid:
                return n
        raise exceptions.MotifGraphException("cannot find node with id")

    #MISC               #######################################################
    def _align_motifs_all_motifs(self):
        start = -1
        for k,v in self.aligned.iteritems():
            if v == 0:
                start = k

        if start == -1:
            raise ValueError("cannot find a place to start in rebuilding motif_graph"
                             " from string")

        i = 0
        for n in graph.transverse_graph(self.graph, start, directed=0):

            if n.index == start:
                self.merger.update_motif(n.data)
                continue
            if n.connections[0] is None:
                continue


            i += 1
            c = n.connections[0]
            parent = c.partner(n.index)
            parent_end_index = c.end_index(parent.index)

            m_added = motif.get_aligned_motif(parent.data.ends[parent_end_index],
                                              n.data.ends[0],
                                              n.data)
            n.data = m_added
            self.merger.update_motif(n.data)

    def _steric_clash(self, m):
        beads = m.beads
        for n in self.graph:
            for c1 in n.data.beads:
                for c2 in beads:
                    if c1.btype == residue.BeadType.PHOS or \
                                    c2.btype == residue.BeadType.PHOS:
                        continue
                    dist = util.distance(c1.center, c2.center)
                    if dist < self.clash_radius:
                        return 1
        return 0

    def _has_ideal_helices(self):
        for n in self.graph:
            if n.data.mtype == motif_type.HELIX and len(n.data.residues()) > 4:
                return 1
        return 0




























