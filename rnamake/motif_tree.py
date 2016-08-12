import base
import option
import settings
import motif
import util
import residue
import motif_type
import tree
import resource_manager as rm
import secondary_structure_factory as ssfactory
import secondary_structure
import motif_merger
import motif_connection
import exceptions

from collections import defaultdict

def motif_tree_from_topology_str(s):
    mt = MotifTree()
    spl = s.split("|")
    node_spl = spl[0].split()
    for i, e in enumerate(node_spl):
        n_spl = e.split(",")
        m = rm.manager.get_motif(name=n_spl[0], end_name=n_spl[1], end_id=n_spl[2])
        if i == 0:
            mt.add_motif(m)
        else:
            mt.add_motif(m, parent_index=int(n_spl[3]), parent_end_index=int(n_spl[4]))

    connection_spl = spl[1].split()
    for c_str in connection_spl:
        c_spl = c_str.split(",")
        mc = motif_connection.MotifConnection(int(c_spl[0]), int(c_spl[1]),
                                              c_spl[2], c_spl[3])

    return mt


def motif_tree_from_ss_tree(sst):
    mt = MotifTree()
    for n in sst.tree:
        #mt.add_motif()
        if n.data.mtype == motif_type.HELIX:
            mt.add_motif(m_name=n.data.name, parent_index=n.parent_index(),
                         parent_end_index=n.parent_end_index())
        else:
            raise ValueError("not implemented")

    return mt


class MotifTree(base.Base):
    """
    MotifTree class orchestrates the connection of motifs to both other motifs
    and full structure and is a core feature of this package

    :attributes:

    `tree` : Tree object
        Holds the motifs in tree form
    `options` : Options object
        Hold various options that effect the behavior of the MotifTree
    `merger` : MotifTreeMerger object
        merges motifs together into a single motif object

    :examples:

    .. code-block:: python

        >>> mt = MotifTree()
        >>> mlib = MotifLibrary()
        >>> mt.add_motif(mlib.get_motif("HELIX.IDEAL"))
        >>> mt.add_motif(mlib.get_motif("HELIX.IDEAL"))

    """

    class _MotifTreePrinter(object):
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

        def __init__(self, mt):
            self.mt = mt
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
            for n in self.mt:
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
            for n in self.mt:
                self.nodes_per_level[self.levels[n.index]] += 1

            for i, n in enumerate(self.mt):
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
                else:
                    pass

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
                if n.parent is not None:
                    parent_end_index = n.parent_end_index()
                    parent_end_name = n.parent.data.ends[parent_end_index].name()
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
            for n in self.mt:
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


    def __init__(self, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.clash_radius = settings.CLASH_RADIUS

        self.tree = tree.TreeStatic()
        self.merger = motif_merger.MotifMerger()
        self.connections = []

    def __len__(self):
        return len(self.tree)

    def __iter__(self):
        self.tree.__iter__()
        return self

    def __repr__(self):
        s = "<(MotifTree: #nodes: %d\n" % (len(self.tree))
        for n in self.tree:
            c_str = ""
            for c in n.children:

                if c is not None:
                    c_str += str(c.index) + ", "
            s += "\t index: %d name: %s end_name: %s children: %s\n" % (n.index, n.data.name, n.data.ends[0].name(), c_str)
        s += ")"

        return s

    def next(self):
        return self.tree.next()

    def copy(self):
        """
        generates a deep copy of the current motif tree

        :return: deep of copy of current tree
        :rtype: MotifTree
        """

        mt = MotifTree()
        new_tree = self.tree.copy()
        mt.tree = new_tree
        mt.merger = self.merger.copy([n.data for n in new_tree.nodes])
        mt.connections = [c.copy() for c in self.connections]
        return mt

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1}

        self.options = option.Options(options)
        self.constraints = {}

    def _validate_arguments_to_add_motif(self, m, m_name):
        """
        makes sure the add_motif function is called correctly

        :param m: motif to add to tree
        :type m: Motif object

        :param m_name: name of motif to add
        :type m_name: str

        :return: None
        """

        if m is not None and m_name is not None:
            raise exceptions.MotifTreeException(
                "cannot supply both a motif and motif name to add a motif to "
                "a motif tree")

        if m is None and m_name is None:
            raise exceptions.MotifTreeException(
                "must supply a motif object or motif name to add_motif")

        if m is not None:
            for n in self.tree.nodes:
                if n.data.id == m.id:
                    raise exceptions.MotifTreeException(
                        "cannot add motif: " + m.name + " to tree as its uuid is " +
                        "already present in the tree")

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
            raise exceptions.MotifTreeException(
                "cannot add motif to tree, motif cannot be found in resource "
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

        parent = self.tree.last_node

        if parent_index != -1:
            try:
                parent = self.tree.get_node(parent_index)
            except exceptions.TreeIndexException:
                raise exceptions.MotifTreeException(
                    "parent_index supplied: " + str(parent_index) + " does not " +
                    "exist in current motif tree")

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

        if parent_end_index != -1 and parent_end_name != "":
            raise exceptions.MotifTreeException(
                "cannot supply parent_end_index and parent_end_name together")

        elif parent_end_name is not None:
            parent_ends = parent.data.get_basepair(name=parent_end_name)
            if len(parent_ends) == 0:
                raise exceptions.MotifTreeException(
                    "cannot find parent_end_name: " + parent_end_name + " in "
                    "parent motif: " + parent.data.name)
            if len(parent_ends) > 1:
                raise exceptions.MotifTreeException(
                    "more then one end was found with parent_end_name: " +
                    parent_end_name + " in parent motif: " + parent.data.name)

            parent_end = parent_ends[0]
            parent_end_index = parent.data.ends.index(parent_end)

            if parent_end_index == parent.data.block_end_add:
                raise exceptions.MotifTreeException(
                    "cannot add motif: to tree as the parent_end_name" +
                    " supplied is blocked see class Motif")

            available = parent.available_pos(parent_end_index)
            if not available:
                raise exceptions.MotifTreeException(
                    "cannot add motif to tree as the end " +
                    "you are trying to add it to is already filled or does "
                    "not exist")

            return [parent_end_index]

        elif parent_end_index != -1:
            if parent_end_index == parent.data.block_end_add:
                raise exceptions.MotifTreeException(
                    "cannot add motif: to tree as the parent_end_index" +
                    " supplied is blocked see class Motif")

            available = parent.available_pos(parent_end_index)
            if not available:
                raise exceptions.MotifTreeException(
                    "cannot add motif to tree as the end " +
                    "you are trying to add it to is already filled or does "
                    "not exist")
            return [parent_end_index]

        else:
            avail_pos = parent.available_children_pos()
            avail_pos.remove(0)
            return avail_pos

    def _steric_clash(self, m):
        for n in self.tree:
            result = motif.clash_between_motifs(n.data, m)
            if result == 1:
                return 1
        return 0

    def add_motif(self, m=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None, m_name=None, m_end_name=None):

        """
        interface to add a motif to the current tree.

        :param m: a motif object to add to tree
        :type m: Motif object

        :param parent_index: the index associated with the parent you want to
            align to. Default is -1 which will use the last motif added if one
            has been added.
        :type parent_index: int

        :param parent_end_index: the end position to align the current motif to
            on the parent motif. Default is -1 or not specified will use first
            end available
        :type parent_end_index: int

        :param parent_end_name:

        """

        self._validate_arguments_to_add_motif(m, m_name)

        parent = self._get_parent_node(parent_index)

        if m is None and m_name is not None:
            m = self._get_motif_from_manager(m_name, m_end_name)

        if parent is None:
            m_copy = m.copy()
            m_copy.get_beads([m_copy.ends[0]])
            self.merger.add_motif(m_copy)
            return self.tree.add_data(m_copy, len(m_copy.ends), -1, -1)

        avail_pos = self._get_parent_available_ends(parent, parent_end_index,
                                                    parent_end_name)

        for p in avail_pos:
            m_added = motif.get_aligned_motif(parent.data.ends[p], m.ends[0], m)
            if self.option('sterics'):
                if self._steric_clash(m_added):
                    continue

            pos =  self.tree.add_data(m_added, len(m_added.ends), parent.index, p)
            self.merger.add_motif(m_added, m_added.ends[0],
                                  parent.data, parent.data.ends[p])
            return pos

        return -1

    def add_motif_tree(self, mt, parent_index=-1, parent_end_name=None):
        parent = self._get_parent_node(parent_index)
        parent_avail_ends = self._get_parent_available_ends(
                                parent, -1, parent_end_name)

        if len(parent_avail_ends) > 1 or len(parent_avail_ends) == 0:
            parent_end_index = -1
        else:
            parent_end_index = parent_avail_ends[0]

        index_hash = {}
        for i, n in enumerate(mt):
            m = rm.manager.get_motif(name=n.data.name, end_name=n.data.ends[0].name())
            if i == 0:
                j = self.add_motif(m, parent_index, parent_end_index)
            else:
                pi = index_hash[n.parent_index()]
                j = self.add_motif(m, pi, parent_end_index=n.parent_end_index())
                if j == -1:
                    self.write_pdbs()
                    print m.name, self.option('sterics')
                    raise ValueError("cannot add_motif_tree")

            index_hash[n.index] = j

    def replace_motif(self, pos, new_motif):
        node = self.get_node(pos)
        if len(new_motif.ends) != len(node.data.ends):
            raise ValueError("attmped to replace a motif with a different number of ends")

        new_motif = new_motif.copy()
        self.merger.replace_motif(node.data, new_motif)
        node.data = new_motif

        for n in tree.transverse_tree(self.tree, pos):
            parent = n.parent
            if parent is None:
                continue
            pei = n.parent_end_index()
            m_added = motif.get_aligned_motif(parent.data.ends[pei],
                                              n.data.ends[0],
                                              n.data)
            n.data = m_added
            self.merger.update_motif(n.data)

    def get_node(self, i):
        return self.tree.get_node(i)

    def write_pdbs(self,name="node"):
        for n in self.tree:
            n.data.to_pdb(name+"."+str(n.index)+".pdb")

    def remove_node(self, i=-1):
        """
        remove a node from the current motif tree. Note there is no checks
        to see whether the node you are removing is a leaf. Thus you could
        be removing a section of the tree. call node.is_leaf() if you would
        like to check this
        """
        if i == -1:
            i = self.tree.last_node.index

        n = self.get_node(i)
        self.tree.remove_node(index=i)
        self.merger.remove_motif(n.data)

    def remove_node_level(self, level=None):
        if level is None:
            level = self.level

        r = range(1, len(self.nodes))
        for i in r[::-1]:
            if self.nodes[i].level >= level:
                self.remove_node(self.nodes[i])
        self.last_node = self.nodes[-1]

    def last_node(self):
        return self.tree.last_node

    def to_pose(self, chain_closure=0):

        pose = self.merger.merge(self.tree)
        self.merger.reset()
        return pose

    def secondary_structure(self):
        return self.merger.secondary_structure()

    def designable_secondary_structure(self):
        ss = self.merger.secondary_structure()

        for n in self.tree.nodes:
            if n.data.name != "HELIX.IDEAL":
                continue
            for r in n.data.residues():
                r_ss = ss.get_residue(uuid=r.uuid)
                if r_ss is not None:
                    r_ss.name = "N"

        return ss

    def to_pdb(self, fname="mt.pdb", renumber=-1, close_chain=0):
        self.merger.get_structure().to_pdb(fname, renumber=renumber,
                                           close_chain=close_chain)

    def to_pdb_str(self, renumber=-1, close_chain=0):
        return self.merger.to_pdb_str(renumber=renumber, close_chain=close_chain)

    def to_pretty_str(self):
        printer = self._MotifTreePrinter(self)
        return printer.print_tree()

    def topology_to_str(self):
        s = ""
        for n in self.tree.nodes:
            s += n.data.name + "," + n.data.ends[0].name() + "," +\
                 n.data.end_ids[0] + "," + str(n.parent_index()) + \
                 "," + str(n.parent_end_index())  +  " "
        s += "|"
        for c in self.connections:
            s += c.to_str() + " "
        return s

    def add_connection(self, i, j, i_bp_name="", j_bp_name=""):
        node_i = self.get_node(i)
        node_j = self.get_node(j)

        node_i_indexes = []
        node_j_indexes = []
        if i_bp_name != "":
            ei = node_i.data.get_end_index(i_bp_name)
            if not node_i.available_pos(ei):
                raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                                 "using bp: " + i_bp_name + "as its not available")
            node_i_indexes.append(ei)
        else:
            node_i_indexes = node_i.available_children_pos()
            node_i_indexes.remove(0)

        if j_bp_name != "":
            ei = node_j.data.get_end_index(j_bp_name)
            if not node_j.available_pos(ei):
                raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                                 "using bp: " + j_bp_name + "as its not available")
            node_j_indexes.append(ei)
        else:
            node_j_indexes = node_j.available_children_pos()
            node_j_indexes.remove(0)

        if len(node_i_indexes) > 1 or len(node_j_indexes) > 1:
            raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                             "its unclear which ends to attach")
        if len(node_i_indexes) == 0 or len(node_j_indexes) == 0:
            raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                             " one node has no available ends")

        self.connections.append(motif_connection.MotifConnection(i, j,
                                                                 i_bp_name, j_bp_name))
        self.merger.connect_motifs(node_i.data, node_j.data,
                                   node_i.data.ends[node_i_indexes[0]],
                                   node_j.data.ends[node_j_indexes[0]])

    def leafs_and_ends(self):
        leaf_nodes = []
        for n in self.tree.nodes:
            f_conn = 0
            for i, c in enumerate(n.children):
                if i == 0:
                    continue
                if c is not None:
                    continue

                leaf_nodes.append([n, i])
        return leaf_nodes

    def residues(self):
        return self.merger.get_structure().residues()

    def get_node_by_id(self, uuid):
        for n in self.tree.nodes:
            if n.data.id == uuid:
                return n
        raise exceptions.MotifTreeException("cannot find node with id")




