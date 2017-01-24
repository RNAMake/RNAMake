import base
import settings
import motif
import util
import residue
import motif_type
import graph
import motif_merger
import copy
import exceptions

from collections import defaultdict

class MotifGraph(object):
    __slots__ = [
        "_clash_radius",
        "_sterics",
        "_graph",
        "_update_merger",
        "_merger",
        "_rm",
        "_aligned",
        "_align_list",
        "_update_align_list"
    ]

    class _MotifGraphBuildPoint(object):
        def __init__(self, node, end_index):
            self.node = node
            self.end_index = end_index

    #SETUP FUNCTIONS ##########################################################
    def __init__(self, rm, mg_str=""):
        self._rm = rm
        self._sterics = 1
        self._graph = graph.GraphStatic()
        self._clash_radius = settings.CLASH_RADIUS
        self._aligned = {}
        self._align_list = []
        self._update_align_list = 1

        self._merger = None
        self._update_merger = 1

        if mg_str != "":
            self._setup_from_str(mg_str)

    def __len__(self):
        return len(self._graph)

    def __iter__(self):
        self._graph.__iter__()
        return self

    def next(self):
        return self._graph.next()

    @classmethod
    def copy(cls, mg):
        new_mg = cls(mg._rm)
        new_mg._graph = mg._graph.copy()
        new_mg._aligned = copy.deepcopy(mg._aligned)
        new_mg._update_merger = 1
        new_mg._update_align_list = 1
        return new_mg

    def _setup_from_str(self, s):
        self._sterics = 0
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

            try:
                rm.manager.get_motif(name=m.name,
                                     end_name=m.ends[0].name())
            except:
                rm.manager.register_motif(m)


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

    def update_indexes(self, index_hash):
        # super hacky, needs to be refactored

        invert_hash = {}
        for new_i, old_i in index_hash.iteritems():
            invert_hash[old_i] = new_i

        largest = 0
        new_aligned = {}
        nodes = []
        for i in invert_hash.iterkeys():
            n = self.get_node(i)
            nodes.append(n)

        for n in nodes:
            new_i = invert_hash[n.index]
            aligned = self._aligned[n.index]
            if new_i > largest:
                largest = new_i

            n.index = new_i
            new_aligned[new_i] = aligned

        self._aligned = new_aligned
        self._update_merger = 1
        self._update_align_list = 1
        self._graph.index = largest+1

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
            for n in self._graph.nodes:
                if n.data.uuid == m.uuid:
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
                m = self._rm.get_motif(name=m_name, end_name=m_end_name)
            else:
                m = self._rm.get_motif(name=m_name)
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

        parent = self._graph.last_node

        if parent_index != -1:
            try:
                parent = self._graph.get_node(parent_index)
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
            parent_end = parent.data.get_basepair(name=parent_end_name)
            if parent_end is None:
                raise exceptions.MotifGraphException(
                    "cannot find parent_end_name: " + parent_end_name + " in "
                    "parent motif: " + parent.data.name)

            parent_end_index = parent.data.get_end_index(parent_end.name)

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
        if parent is None:
            pos = self._graph.add_data(m, -1, -1, -1, m.num_ends(), orphan=1)
            if pos != -1:
                self._aligned[pos] = 0
                self._update_align_list = 1
                self._update_merger = 1
            return pos

        else:
            pos = self._graph.add_data(m, parent.index, parent_end_index,
                                       0, m.num_ends())

            if pos != -1:
                self._aligned[pos] = 1
                self._update_align_list = 1
                self._update_merger = 1
            return pos

    def add_motif(self, m=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None, m_name=None, m_end_name=None,
                  orphan=0):

        self._validate_arguments_to_add_motif(m, m_name)
        parent = self._get_parent_node(parent_index)

        if m is None and m_name is not None:
            m = self._get_motif_from_manager(m_name, m_end_name)
        else:
            if not self._rm.contains_motif(name=m.name,
                                           end_name=m.get_end(0).name):
                self._rm.add_motif(m, m.mtype, m.name)

        if parent is None or orphan:
            return self._add_motif_to_graph(motif.Motif.copy(m, new_uuid=1), None, None)

        avail_pos = self._get_parent_available_ends(parent, parent_end_index,
                                                    parent_end_name)

        for p in avail_pos:
            m_added = motif.get_aligned_motif(parent.data.get_end(p), m.get_end(0), m)
            if self._sterics:
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

        self._graph.connect(i, j, node_i_ei, node_j_ei)
        self._update_merger = 1

    #REMOVE FUNCTIONS   #######################################################
    def remove_motif(self, pos=-1):
        #TODO must update alignments! if something is being aligned and its parent
        # is possible it is no longer being aligned!
        if pos == -1:
            pos = self.last_node().index
        n = self._graph.get_node(pos)
        self._graph.remove_node(pos)
        del self._aligned[pos]
        self._update_merger = 1
        self._update_align_list = 1

    def remove_node_level(self, level=None):
        if level is None:
            level = self._graph.level

        r = range(1, len(self._graph.nodes))
        for i in r[::-1]:
            if self._graph.nodes[i].level >= level:
                self.remove_motif(self._graph.nodes[i].index)

    #DESIGNING          #######################################################
    def designable_secondary_structure(self):
        self.__update_merger()
        ss = self._merger.get_merged_secondary_structure()

        for n in self._graph.nodes:
            if n.data.name != "HELIX.IDEAL":
                continue
            for r in n.data.iter_res():
                r_ss = ss.get_residue(uuid=r.uuid)
                if r_ss is not None:
                    r_ss.set_name("N")

        return ss

    def replace_ideal_helices(self):
        found = 1
        while found:
            found = 0
            for n in self._graph.nodes:
                if n.data.mtype != motif_type.HELIX:
                    continue
                if n.data.num_res() == 4:
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
                old_n_aligned = self._aligned[n.index]
                self.remove_motif(i)
                old_n = n

                h = self._rm.get_motif(name="HELIX.IDEAL")
                if parent is None:
                    pos = self._add_motif_to_graph(h, None, None)

                    if old_n_aligned == 0:
                        # fix alignment for parentless nodes
                        n = self.get_node(pos)
                        m_added = motif.get_aligned_motif(old_n.data.get_end(0),
                                                          n.data.get_end(0),
                                                          n.data)
                        n.data = m_added


                else:
                    m_added = motif.get_aligned_motif(parent.data.get_end(parent_end_index),
                                                      h.get_end(0), h)
                    pos = self._add_motif_to_graph(m_added, parent, parent_end_index)


                for j in range(0, count):
                    h = self._rm.get_motif(name="HELIX.IDEAL")
                    parent =  self._graph.get_node(pos)
                    m_added = motif.get_aligned_motif(parent.data.get_end(1), h.get_end(0), h)
                    pos = self._add_motif_to_graph(m_added, parent, 1)

                if other:
                    self._graph.connect(pos, other.index, 1, other_end_index)
                    self._aligned[other.index] = 1

                break

        self._update_merger = 1
        self._update_align_list = 1
        self._align_motifs_all_motifs()

    def replace_helix_sequence(self, ss=None, seq=None):
        if ss is None:
            ss = self.designable_secondary_structure()
            ss.replace_sequence(seq)

        for n in self._graph.nodes:
            if n.data.mtype != motif_type.HELIX:
                continue
            ss_m = ss.get_motif(n.data.uuid)
            if ss_m.get_end_id(0) == n.data.get_end_id(0) and n.data.name != "HELIX.IDEAL":
                continue

            m = self._rm.get_bp_step(ss_m.get_end_id(0), n.data)

            if self._aligned[n.index] == 0:
                m_added = motif.get_aligned_motif(n.data.get_end(0), m.get_end(0), m)
                n.data = m_added
            else:
                n.data = m

        self._update_merger = 1
        self._align_motifs_all_motifs()

    def replace_motif(self, pos, new_motif):
        node = self.get_node(pos)
        if len(new_motif.ends) != len(node.data.ends):
            raise exceptions.MotifGraphException(
                "attempted to replace a motif with a different number of ends")

        node.data = new_motif.copy()

        self._update_merger = 1
        self._align_motifs_all_motifs()

    #GRAPH WRAPPER      #######################################################
    def increase_level(self):
        self._graph.increase_level()

    def decrease_level(self):
        self._graph.decrease_level()

    def last_node(self):
        return self._graph.last_node

    def get_node(self, i=None, uuid=None, m_name=None):
        if i is not None:
            return self._graph.get_node(i)

        node = None
        for n in self._graph.nodes:
            if n.data.uuid== uuid:
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
        self.__update_merger()
        return self._merger.get_merged_secondary_structure()

    def get_structure(self):
        self.__update_merger()
        return self._merger.get_merged_structure()

    def sequence(self):
        self.__update_merger()
        return self._merger.get_merged_structure().sequence()

    def dot_bracket(self):
        self.__update_merger()
        return self._merger.get_merged_structure().dot_bracket()

    #OUTPUTING          #######################################################
    def nodes_to_pdbs(self, name="node"):
        for n in self._graph.nodes:
            n.data.to_pdb(name + "." + str(n.index) + ".pdb")

    def to_pdb(self, name="test.pdb", renumber=-1, close_chain=0):
        self.__update_merger()
        return self._merger.get_merged_structure().to_pdb(
                    name, renumber=renumber, close_chain=close_chain)

    def topology_to_str(self):
        s = ""
        con_str = ""
        seen_connections = {}
        for n in self._graph.nodes:
            s += n.data.name + "," + n.data.get_end(0).name + "," + str(n.index) + ","
            s += str(self._aligned[n.index]) + "|"
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
        for n in self._graph.nodes:
            s += n.data.to_str() + "^" + str(n.index) + "^"
            s += str(self._aligned[n.index]) + " KAK "
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
        printer = MotifGraphPrinter(self)
        return printer.print_graph()

    #GETTERS            #######################################################
    def get_build_points(self):
        build_points = []
        for n in self._graph.nodes:
            for i, c in enumerate(n.connections):
                if c is None and i != 0:
                    build_points.append(self._MotifGraphBuildPoint(n, i))

        return build_points

    def get_end(self, pos=-1, m_name="", m_end_name=""):
        n = None
        if pos != -1:
            n = self._graph.get_node(pos)
        elif m_name != "":
            nodes = []
            for n in self._graph.nodes:
                if n.data.name == m_name:
                    nodes.append(n)
            if len(nodes) > 1:
                raise ValueError("cannot get end, too many motifs match name given "
                                 + m_name)
            n = nodes[0]

        if m_end_name != "":
            bp = n.data.get_basepair(name=m_end_name)
            if bp is None:
                raise ValueError("found motif but " + m_end_name + "is not is a "
                                                                   "basepair in it")
            end_index = n.data.get_end_index(name=m_end_name)
            end = n.data.get_end(end_index)
            # check to see if the position is available
            self._graph.get_availiable_pos(n, end_index)
            return end

        else:
            avail_pos = n.available_children_pos()
            if len(avail_pos) > 1:
                raise ValueError("too many free ends to pick one in get_end")
            if len(avail_pos) == 0:
                raise ValueError("no ends available in get_end")

            end = n.data.get_end(avail_pos[0])
            return end

    def get_not_aligned_nodes(self):
        not_aligned = []
        for n in self._graph:
            if self._aligned[n.index] == 0:
                not_aligned.append(n)
        return not_aligned

    #MISC               #######################################################
    def _get_align_list(self):
        if not self._update_align_list:
            return self._align_list

        non_aligned_nodes = self.get_not_aligned_nodes()
        self._align_list = []
        used_nodes = {}
        for start in non_aligned_nodes:
            open = [start]
            seen_nodes = {}

            while open:
                n = open.pop(0)
                seen_nodes[n] = 1
                if n.index == start.index:
                    self._align_list.append(n)
                    used_nodes[n] = 1
                else:
                    if n.connections[0] is None:
                        continue
                    c = n.connections[0]
                    parent = c.partner(n.index)

                    if parent not in used_nodes:
                        continue

                    self._align_list.append(n)

                used_nodes[n] = 1
                for i, c in enumerate(n.connections):
                    if i == n.data.block_end_add or c is None:
                        continue

                    partner_n = c.partner(n.index)
                    if partner_n in seen_nodes or partner_n in used_nodes:
                        continue

                    # if something goes wrong check this!
                    if c.end_index(partner_n.index) == partner_n.data.block_end_add:
                        open.append(partner_n)
                    elif partner_n.data.num_ends() == 1:
                        open.append(partner_n)


        self._update_align_list = 0
        return self._align_list

    def _align_motifs_all_motifs(self):
        non_aligned_nodes = self.get_not_aligned_nodes()
        align_list = self._align_list

        for n in align_list:
            if n in non_aligned_nodes:
                continue

            parent = n.connections[0].partner(n.index)
            pei = n.connections[0].end_index(parent.index)
            m_added = motif.get_aligned_motif(parent.data.get_end(pei),
                                              n.data.get_end(0),
                                              n.data)
            n.data = m_added

    def _steric_clash(self, m):
        for n in self._graph.nodes:
            result = motif.clash_between_motifs(n.data, m)
            if result == 1:
                return 1
        return 0

    def _has_ideal_helices(self):
        for n in self._graph:
            if n.data.mtype == motif_type.HELIX and len(n.data.residues()) > 4:
                return 1
        return 0

    def __update_merger(self):
        if not self._update_merger:
            return

        self._merger = motif_merger.MotifMerger(self._rm.motif_factory)

        non_aligned_nodes = self.get_not_aligned_nodes()
        align_list = self._get_align_list()
        seen_connections = {}

        for n in align_list:
            if n in non_aligned_nodes:
                self._merger.add_motif(n.data)
                continue

            c = n.connections[0]
            parent = c.partner(n.index)
            parent_end_index = c.end_index(parent.index)
            seen_connections[c] = 1
            self._merger.add_motif(n.data, n.data.get_end(0),
                                  parent.data, parent.data.get_end(parent_end_index))


        for n in align_list:
            for c in n.connections:
                if c is None:
                    continue
                if c in seen_connections:
                    continue
                partner = c.partner(n.index)
                end1 = n.data.get_end(c.end_index(n.index))
                end2 = partner.data.get_end(c.end_index(partner.index))
                self._merger.connect_motifs(n.data, partner.data, end1, end2)
                seen_connections[c] = 1

        self._update_merger = 0


class MotifGraphPrinter(object):
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



















