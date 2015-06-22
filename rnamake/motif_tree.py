import base
import option
import settings
import motif
import util
import residue
import motif_type
import motif_tree_merger
import graph

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
        self.clash_radius = settings.CLASH_RADIUS - 0.1

        self.graph = graph.GraphStatic()
        if m is not None:
            self.add_motif(m)

        self.merger = motif_tree_merger.MotifTreeMerger()

    def __len__(self):
        return len(self.graph)

    def __iter__(self):
        self.graph.__iter__()
        return self

    def next(self):
        return self.graph.next()

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1,
                    'full_beads_first_res' : 1}

        self.options = option.Options(options)
        self.constraints = {}

    def add_motif(self, m=None, parent_index=-1, parent_end_index=-1):
        parent = self.graph.last_node
        if parent_index != -1:
            parent = self.graph.get_node(parent_index)

        if parent is None:
            m_copy = m.copy()
            m_copy.get_beads(m_copy.ends)
            return self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends))

        avail_pos = self.graph.get_availiable_pos(parent, parent_end_index)

        for p in avail_pos:
            m_added = motif.get_aligned_motif(parent.data.ends[p], m.ends[0], m)
            if self.option('sterics'):
                if self._steric_clash(m_added):
                    continue

            for r in m_added.residues():
                r.new_uuid()

            return self.graph.add_data(m_added, parent.index, p, 0, len(m_added.ends))
        #self._update_beads(parent, new_node)

        return -1

    def get_node(self, i):
        return self.graph.get_node(i)

    def write_pdbs(self,name="node"):
        for n in self.graph:
            n.data.to_pdb(name+"."+str(n.index)+".pdb")

    def remove_node(self, i=-1):
        """
        remove a node from the current motif tree. Note there is no checks
        to see whether the node you are removing is a leaf. Thus you could
        be removing a section of the tree. call node.is_leaf() if you would
        like to check this
        """
        if i == -1:
            i = len(self.graph)-1

        self.graph.remove_node(i)

    def remove_node_level(self, level=None):
        if level is None:
            level = self.level

        r = range(1, len(self.nodes))
        for i in r[::-1]:
            if self.nodes[i].level >= level:
                self.remove_node(self.nodes[i])
        self.last_node = self.nodes[-1]

    def leafs(self):
        leafs = []
        for n in self.graph:
            if len(n.available_children_pos()) > 0:
                leafs.append(n)
        return leafs

    def to_pose(self, chain_closure=0):

        pose = self.merger.merge(self.graph)
        self.merger.reset()
        return pose

    def to_pdb(self, fname="mt.pdb", include_head=1, chain_closure=1):
        pose = self.to_pose()
        pose.to_pdb(fname)

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
        for end in parent.motif.ends:
            if end not in not_used:
                exclude.append(end)

        parent.motif.get_beads(exclude)
        child.motif.get_beads()

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

    def _find_other_connections_to_head(self):
        leafs = self.leafs()
        for leaf in leafs:
            result = self._add_connection(self.nodes[0], leaf)
            if not result:
                continue
            head_node_open_ends = self.nodes[0].get_available_ends()
            if len(head_node_open_ends) == 0:
                break

    def _add_connection(self, node_1, node_2, cutoff=25):
        if node_1 == node_2:
            return 0

        avail_ends_1 = node_1.available_ends()
        avail_ends_2 = node_2.available_ends()

        for end1 in avail_ends_1:
            for end2 in avail_ends_2:
                dist = util.distance(end1.d(), end2.d())
                if dist < cutoff:
                    new_connection = MotifTreeConnection(node_1, node_2, end1,
                                                         end2)







