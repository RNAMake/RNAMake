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
    MotifTree class orchestrates the connection of motifs to both other motifs and full structure and is a core feature of this package

    :param motif: the motif that will serve as the head node that everything will be built from, if nothing is specified a single basepair will be used
    :type motif: Motif object

    Examples
m
    .. code-block:: python
        >>>mt = MotifTree()
        >>>mlib = MotifLibrary()
        >>>mt.add_motif(mlib.get_motif("HELIX.IDEAL"))
        >>>mt.add_motif(mlib.get_motif("HELIX.IDEAL"))

        #number of nodes added
        >>>mt.nodes
        3

        #merge motifs together and get a single motif and you can print that a pdb
        >>>mm = mt.get_merged_motif()
        >>>mm.to_pdb("test.pdb")

    Attributes
    ----------
    `nodes` : List of MotifTreeNodes object
        All the nodes that belong to this MotifTree, each node contains the connection information between nodes and the motif
    `clash_radius` : Float
        the radius between beads that will stop a motif from being added
    `level` : Int
        The level of nodes being built, this permits for fast removal of bulk nodes
    `options` : Options object
        Hold various options that effect the behavior of the MotifTree
    `last_node` : MotifTreeNode
        The last node added to the tree useful for quick additions to the tree
    `merger` : MotifTreeMerger object
        merges motifs together into a single motif object

    """

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

    def add_motif(self, m=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None, m_name=None, m_end_name=None):

        if m is None and m_name is not None:
            if m_end_name is not None:
                m = rm.manager.get_motif(name=m_name, end_name=m_end_name)
            else:
                m = rm.manager.get_motif(name=m_name)

        for n in self.tree.nodes:
            if n.data.id == m.id:
                raise exceptions.MotifTreeException(
                    "cannot add motif: " + m.name + " to graph as its uuid is "
                    "already present in the graph")

        parent = self.tree.last_node
        if parent_index != -1:
            parent = self.tree.get_node(parent_index)

        if parent is None:
            m_copy = m.copy()
            m_copy.get_beads(m_copy.ends)
            self.merger.add_motif(m_copy)
            return self.tree.add_data(m_copy, len(m_copy.ends), -1, -1)

        if parent_end_name is not None:
            parent_end = parent.data.get_basepair(name=parent_end_name)[0]
            parent_end_index = parent.data.ends.index(parent_end)

        avail_pos = self.tree.get_available_pos(parent, parent_end_index)

        #print avail_pos
        for p in avail_pos:
            if p == parent.data.block_end_add:
                continue

            m_added = motif.get_aligned_motif(parent.data.ends[p], m.ends[0], m)
            if self.option('sterics'):
                if self._steric_clash(m_added):
                    continue

            pos =  self.tree.add_data(m_added, len(m_added.ends), parent.index, p)
            self.merger.add_motif(m_added, m_added.ends[0],
                                  parent.data, parent.data.ends[p])
            return pos

        return -1

    def add_motif_tree(self, mt, parent_index=-1, parent_end_name=""):
        if parent_index == -1:
            for n in mt:
                self.add_motif(n.data)
            return

        parent = self.get_node(parent_index)
        bps = parent.data.get_basepair(name=parent_end_name)
        if len(bps) == 0:
            raise ValueError("cannot find parent end in add_motif_tree")
        pei = parent.data.ends.index(bps[0])

        for i, n in enumerate(mt):
            if i == 0:
                self.add_motif(n.data, parent_index, pei)
            else:
                self.add_motif(n.data)

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

    def _steric_clash(self, m):
        beads = m.beads
        for n in self.tree:
            for c1 in n.data.beads:
                for c2 in beads:
                    if c1.btype == residue.BeadType.PHOS or \
                       c2.btype == residue.BeadType.PHOS:
                        continue
                    dist = util.distance(c1.center, c2.center)
                    if dist < self.clash_radius:
                        print dist
                        return 1
        return 0

    def _find_other_connections_to_head(self):
        leafs = self.leafs()
        for leaf in leafs:
            result = self._add_connection(self.nodes[0], leaf)
            if not result:
                continue
            head_node_open_ends = self.nodes[0].get_available_ends()
            if len(head_node_open_ends) == 0:
                break

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

        #self.graph.connect(i, j, node_i_indexes[0], node_j_indexes[0])

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




