import directed_graph
from . import exceptions

class MotifType(object):
    ALL_ATOM   = 0
    STATE      = 1
    SSTRUCTURE = 2


class MotifTypeDirectedGraph(object):
    __slots__ = [
        "_dg",
        "_sterics",
        "_motif_class_type",
        "_motif_aligner",
        "_motif_type",
        "_motif_merger",
    ]

    def __init__(self):
        self._dg = directed_graph.DirectedGraph()
        # should we worry about residues colliding
        self._sterics = 1
        # what type of motifs are these
        self._motif_class_type = None
        # how do we align these motifs to each other
        self._motif_aligner = None

    def __len__(self):
        return len(self._dg)

    def __iter__(self):
        return self._dg.__iter__()

    @classmethod
    def copy(cls, mdg):
        new_mdg = cls(mdg._rm)
        # add all motif to new graph
        for i in mdg:
            m = mdg.get_motif(i)
            m_copy = new_mdg._motif_class_type.copy(m)
            if mdg.has_parent(i):
                new_mdg.add_motif(m_copy, mdg.get_parent_index(i),
                                  mdg.get_parent_end_index(i))
            else:
                new_mdg.add_motif(m)

        edges = mdg.get_all_edges()
        for full_edge in edges:
            ni, nj, ei, ej = full_edge
            if not new_mdg.connection_exists(ni, nj, ei, ej):
                new_mdg.add_connection(ni, nj, ei, ej)
        return new_mdg

    def _from_str(cls, s, rts):
        pass

    def get_motif(self, i=None, uuid=None):
        if i is not None:
            return self._dg.get_node(i)

    def _get_parent_end_index(self, parent_index, parent_end_index, parent_end_name):
        if parent_index == -1:
            return -1
        # must supply either an end index or end name
        if parent_end_index != -1 and parent_end_name is not None:
            raise ValueError("cannot specify parent_end_index and parent_end_name")
        # can't supply both
        if parent_end_index == -1 and parent_end_name is None:
            raise ValueError("must supply either parent_end_index or parent_end_name")

        p = self._dg.get_node(parent_index)
        pi = parent_end_index

        if parent_end_name is not None:
            try:
                pi = p.get_end_index(parent_end_name)
            except exceptions.RNAStructureException as e:
                raise ValueError(e)

        if pi == p.block_end_add:
            raise ValueError("cannot add motif at block_end_add position")

        return pi

    def _is_motif_in_graph_already(self, m):
        for ni in self._dg:
            n_motif = self._dg.get_node(ni)
            if m.uuid == n_motif.uuid:
                raise ValueError(
                    "cannot add motif: " + m.name + " to graph as its uuid is " +
                    "already present in the graph")
        return 0

    def add_motif(self, m, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None, index=-1):

        # check to see if this specific motif object already exists in graph
        self._is_motif_in_graph_already(m)
        # handle option arguments and conver parent_end_name to index
        parent_end_index = self._get_parent_end_index(
                                parent_index, parent_end_index, parent_end_name)
        # get aligned motif
        m_aligned = self._align(parent_index, parent_end_index, m)
        # check sterics
        if self._steric_clash(m_aligned):
            raise ValueError(
                "cannot add motif: "+ m.name + " it has steric clashes at "
                "position: " + str(self._dg.index))
        # add motif to graph
        return self._dg.add_node(m_aligned, m.num_ends(), parent_index,
                                 parent_end_index, 0, index)

    def remove_motif(self, i):
        self._dg.remove_node(i)
        self._update_motif_alignments()

    def _align(self, parent_index, parent_end_index, m):
        m_copy = self._motif_class_type.copy(m)
        if parent_index == -1:
            return m_copy
        parent = self._dg.get_node(parent_index)
        parent_end = parent.get_end(parent_end_index)
        if parent_end is None:
            raise ValueError(
                "cannot add motif to graph parent does not have end specified")
        self._motif_aligner.align(parent_end, m_copy)
        return m_copy

    def _steric_clash(self, m):
        if not self._sterics:
            return 0
        for ni in self._dg:
            if m.steric_clash(self._dg.get_node(ni)):
                return 1
        return 0

    def _get_parent_end(self, i):
        pei = self._dg.get_parent_edge_index(i)
        return self._dg.get_parent(i).get_end(pei)

    def _update_motif_alignments(self, start_i=-1):
        start = 0
        if start_i == -1:
            start = 1

        for i in self._dg:
            if not start:
                if start_i == i:
                    start = 1
                    continue
            else:
                if not self._dg.has_parent(i):
                    continue

                self._motif_aligner.align(self._get_parent_end(i),
                                          self._dg.get_node(i))

    def _get_motif_from_manager(self, m_name, m_end_name):
        raise NotImplementedError

    # wrappers from directed graph
    def get_parent_index(self, i):
        return self._dg.get_parent_index(i)

    def get_parent_end_index(self, i):
        return self._dg.get_parent_edge_index(i)

    def has_parent(self, i):
        return self._dg.has_parent(i)

    def get_all_edges(self):
        return self._dg.get_all_edges()

    def are_motif_connected(self, ni, nj):
        return self._dg.are_nodes_connected(ni, nj)

    def get_end_indexes(self, ni, nj):
        return self._dg.get_edge_indexes(ni, nj)

    def add_connection(self, ni, nj, ei, ej):
        return self._dg.add_edge(ni, nj, ei, ej)

    def connection_exists(self, ni, nj, ei, ej):
        return self._dg.edge_exists(ni, nj, ei, ej)

    def is_end_filled(self, ni, ei):
        return self._dg.edge_index_empty(ni, ei)

    @property
    def sterics(self):
        return self._sterics

    @sterics.setter
    def sterics(self, new_value):
        if new_value != 0 and new_value != 1:
            raise ValueError("sterics can only be set to 0 or 1")
        self._sterics = new_value

    @property
    def last(self):
        return self._dg.last

    @property
    def last_motif(self):
        return self.get_motif(self.last)