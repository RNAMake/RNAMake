from . import base
from . import option
from . import settings
from . import motif
from . import util
from . import residue
from . import motif_type
from . import motif_tree_merger

class MotifTree(base.Base):
    """
    MotifTree class orchestrates the connection of motifs to both other motifs and full structure and is a core feature of this package

    :param motif: the motif that will serve as the head node that everything will be built from, if nothing is specified a single basepair will be used
    :type motif: Motif object

    Examples

    .. code-block:: python
        >>>mt = MotifTree()
        >>>mlib = MotifLibrary()
        >>>mt.add_motif(mlib.get_motif("HELIX.IDEAL))
        >>>mt.add_motif(mlib.get_motif("HELIX.IDEAL))

        #number of nodes added
        >>>mt.nodes
        3

        #merge motifs together and get a single motif and you can print that a pdb
        >>>mm = mt.get_merged_motif()
        >>>mm.to_pdb("test.pdb)

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

    def __init__(self, m=None, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)

        if m is None:
            head = self._get_default_head()
        else:
            head = MotifTreeNode(m.copy(), 0, 0, 0)
            if self.option('sterics') == 1:
                head.motif.get_beads(head.motif.ends)

        self.nodes = [ head ]
        self.clash_radius = 2.5
        self.level = 1
        self.last_node = head
        self.merger = motif_tree_merger.MotifTreeMerger()

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1,
                    'full_beads_first_res' : 1,
                    'ideal_bp_score'       : -2 # TODO fix
                    }

        self.options = option.Options(options)
        self.constraints = {}

    def add_motif(self, m, end_index=None, end=None, end_flip=None,
                  parent=None, parent_index=None, parent_end=None):

        # if parent is not specified use the last node added
        if parent is None:
            parent = self.last_node

        # get which ends are available to align
        parent_ends = self._parse_ends(parent.motif, parent_index, parent_end,
                                       parent)
        ends = self._parse_ends(m, end_index, end)

        if len(parent_ends) == 0:
            raise ValueError("cannot add to this parent node, it has no \
                             available ends")

        if len(ends) == 0:
            raise ValueError("cannot add this node to the tree, it has no \
                             available ends")

        new_node = None
        flip_status = [x for x in range(len(m.ends))]
        if end_flip is not None:
            flip_status = [end_flip]
        for parent_end in parent_ends:
            for end in ends:
                if new_node is not None:
                    break
                new_node = self._add_motif(parent, m, parent_end, end,
                                           flip_status)
        if new_node is not None:
            self.nodes.append(new_node)
            self.last_node = new_node
        return new_node

    def write_pdbs(self,name="node"):
        for i,n in enumerate(self.nodes):
            n.motif.to_pdb(name+"."+str(i)+".pdb")

    def remove_node(self, node):
        """
        remove a node from the current motif tree. Note there is no checks
        to see whether the node you are removing is a leaf. Thus you could
        be removing a section of the tree. call node.is_leaf() if you would
        like to check this

        :param node: node object that you would like to remove from the tree
        :type node: MotifTreeNode object

        """
        if node.index == 0:
            raise ValueError("cannot remove head node from tree")

        connected = node.connected_nodes()

        while 0 < len(node.connections):
            node.connections[-1].disconnect()

        self.last_node = connected[0]
        self.nodes.remove(node)

    def leafs(self):
        leafs = []
        for n in self.nodes:
            if n.is_leaf():
                leafs.append(n)
        return leafs

    def get_pose(self, include_head=0, chain_closure=0):
        if include_head:
            self._find_other_connections_to_head()

        return self.merger.merge(self, include_head=include_head,
                                 chain_closure=chain_closure)

    def to_pdb(self, fname="mt.pdb", include_head=0, chain_closure=0):
        pose = self.get_pose(include_head=include_head,
                             chain_closure=chain_closure)
        pose.to_pdb(fname)

    def to_str(self):
        s = ""
        for n in self.nodes:
            s += n.to_str() + "#"
        return s

    def _parse_ends(self, m, end_index, end_bp, node=None):
        ends = []
        if   end_index is not None and end_bp is not None:
            raise ValueError("MotifTree: Cannot supply both end_index and \
                             end_bp to add_motif")
        elif end_index is not None:
            try:
                ends = [ m.ends[end_index] ]
            except:
                raise ValueError("end index: " + end_index + " is out of \
                                 out of bounds in add_motif")
        elif end_bp is not None:
            if end_bp not in m.ends:
                raise ValueError("end_bp is not an end of motif in add_motif")
            ends = [ end_bp ]
        else:
            if node:
                ends = node.available_ends()
            else:
                ends = m.ends

        return ends

    def _add_motif(self, parent, m, parent_end, end, flip_status=[0,1]):
        for flip in flip_status:
            end.flip(flip)
            motif.align_motif(parent_end, end, m)

            if self.option('sterics') == 1:
                if self._steric_clash(m, end):
                    m.reset()
                    continue

            new_node = MotifTreeNode(m.copy(), self.level, len(self.nodes),
                                     flip)

            m.reset()
            new_end = new_node.motif.get_basepair(uuid1=end.res1.uuid,
                                                  uuid2=end.res2.uuid)[0]
            MotifTreeConnection(parent, new_node, parent_end, new_end)
            #if new uuids are not assigned there will be duplicate
            #uuid after adding multiple bp steps
            for res in new_node.motif.residues():
                res.new_uuid()

            if len(self.nodes) == 1 and self.option('full_beads_first_res') and\
               self.option('sterics'):
                new_node.motif.get_beads()

            self._update_beads(parent, new_node)
            return new_node
        return None

    def _update_beads(self, parent, child):
        """
        This may seem strange but it correct the small differences that can
        occur between the beads in a motif tree vs that in a merged motif.
        Merged motifs remove overlapping basepairs by including the motif's
        basepair over the helices basepair. Thus this function reflects this
        change a motif tree will behave like a merged motif sterically.
        """

        if parent.motif.mtype != motif_type.HELIX or \
           child.motif.mtype != motif_type.HELIX:
            return

        not_used = parent.available_ends()
        exclude = []
        for end in child.motif.ends:
            if end not in not_used:
                exclude.append(end)

        parent.motif.get_beads(exclude)
        child.motif.get_beads()

    def _get_default_head(self):
        mdir = settings.RESOURCES_PATH + "/start"
        m = motif.Motif(mdir)
        return MotifTreeNode(m, 0, 0, 0)

    def _steric_clash(self, m, end):
        beads = m.get_beads([end])

        for n in self.nodes[::-1]:
            for c1 in n.motif.beads:
                for c2 in beads:
                    if c1.btype == residue.BeadType.PHOS or \
                       c2.btype == residue.BeadType.PHOS:
                        continue
                    dist = util.distance(c1.center, c2.center)
                    if dist < self.clash_radius:
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

    def _add_connection(self, node_1, node_2):
        if node_1 == node_2:
            return 0

        avail_ends_1 = node_1.available_ends()
        avail_ends_2 = node_2.available_ends()

        for end1 in avail_ends_1:
            for end2 in avail_ends_2:
                dist = util.distance(end1.d(), end2.d())
                if dist < 10:
                    new_connection = MotifTreeConnection(node_1, node_2, end1,
                                                         end2)


class MotifTreeNode(object):
    def __init__(self, m, level, index, flip):
        self.motif, self.level, self.index, self.flip = m, level, index, flip
        self.end_status = { end.uuid : 1 for end in self.motif.ends }
        self.connections = []

    def available_ends(self):
        ends = []
        for uuid in self.end_status.iterkeys():
            for end in self.motif.ends:
                if end.uuid == uuid:
                    ends.append(end)
                    break
        ends.sort(key=lambda x : self.motif.ends.index(x), reverse=True)
        return ends

    def connected_nodes(self):
        return [c.partner(self) for c in self.connections]

    def is_leaf(self):

        #the head is not a leaf
        if self.index == 0:
            return 0

        if len(self.connected_nodes()) == 1:
            return 1
        else:
            return 0

    def parent(self):
        for n in self.connected_nodes():
            if n.index < self.index:
                return n

    def connection(self, n):
        for c in self.connections:
            if c.node_1 == n or c.node_2 == n:
                return c

    def to_str(self):
        s = self.motif.to_str() + "!"
        parent = self.parent()
        if parent is None:
            s += "-1!-1!-1!-1"
            return s

        c = self.connection(parent)
        parent_end = c.motif_end(parent)
        node_end = c.motif_end(self)
        s += str(parent.index) + "!" + str(parent.motif.ends.index(parent_end)) +\
             "!" + str(self.motif.ends.index(node_end)) + "!" + str(self.flip)
        return s


class MotifTreeConnection(object):
    """
    A connection between two MotifTreeNodes in MotifTree class. Connections
    represent basepair alignment overlap.

    :param node_1: First node to be connected
    :param node_2: Second node to be connected
    :param end_1: Basepair end from node_1 that is aligned with the end_2 from
        node_2
    :param end_2: Basepair end from node_2 that is aligned with the end_1 from
        node_1
    :param no_overlap: controls whether both end_1 and end_2 should be kept
        instead of keeping of only one

    :type node_1: MotifTreeNode object
    :type node_2: MotifTreeNode object
    :type end_1: Basepair object
    :type end_2: Basepair object
    :type no_overlap: Int

    Attributes
    ----------
    `node_1` : MotifTreeNode object
        First node to be connected
    `node_2` : MotifTreeNode object
        Second node to be connected
    `end_1` : Basepair object
        Basepair end from node_1 that is aligned with the end_2 from node_2
    `end_2` : Basepair object
        Basepair end from node_2 that is aligned with the end_1 from node_1
    `no_overlap` : Bool
        Controls whether both end_1 and end_2 should be kept instead of keeping
        of only one
    """

    def __init__(self, node_1, node_2, end_1, end_2, no_overlap=0):
        self.node_1, self.node_2, self.end_1, self.end_2, self.no_overlap = \
            node_1, node_2, end_1.uuid, end_2.uuid, no_overlap

        self.node_1.end_status[end_1] = 0
        self.node_2.end_status[end_2] = 0
        self.node_1.connections.append(self)
        self.node_2.connections.append(self)

    def disconnect(self):
        self.node_1.end_status[self.end_1] = 1
        self.node_2.end_status[self.end_2] = 1
        self.node_1.connections.remove(self)
        self.node_2.connections.remove(self)

    def partner(self, node):
        if   node == self.node_1:
            return self.node_2
        elif node == self.node_2:
            return self.node_1
        else:
            raise ValueError("node is not in connection object cannot call \
                             pair")

    def motif_end(self, node):
        if   node == self.node_1:
            return self.node_1.motif.get_basepair(self.end_1)[0]
        elif node == self.node_2:
            return self.node_2.motif.get_basepair(self.end_2)[0]
        else:
            raise ValueError("node is not in connection object cannot call \
                             motif_end")


def str_to_motif_tree(s):
    spl = s.split("#")
    for i, e in enumerate(spl[:-1]):
        mtn_spl = e.split("!")
        m = motif.str_to_motif(mtn_spl[0])
        if i == 0:
            mt = MotifTree(m)
            continue
        mtn = MotifTreeNode(m, mt.level, int(mtn_spl[3]), int(mtn_spl[4]))
        mtn.index = i
        parent = mt.nodes[ int(mtn_spl[1]) ]
        parent_end = parent.motif.ends[ int(mtn_spl[2]) ]
        node_end = m.ends[ int(mtn_spl[3]) ]
        MotifTreeConnection(parent, mtn, parent_end, node_end)
        mt.nodes.append(mtn)
    return mt







