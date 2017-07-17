import directed_graph
import motif

class MotifTypeDirectedGraph(object):
    __slots__ = [
        "_dg",
        "_sterics",
        "_motif_class_type"
    ]

    def __init__(self):
        self._dg = directed_graph.DirectedGraph()
        self._sterics = 1
        self._motif_class_type = motif.Motif

    def __len__(self):
        return len(self._dg)

    def __iter__(self):
        return self._dg.__iter__()

    def _get_parent_end_index(self, parent_index, parent_end_index, parent_end_name):
        if parent_index == -1:
            return -1

        if parent_end_index != -1 and parent_end_name is not None:
            raise ValueError("cannot specify parent_end_index and parent_end_name")

        if parent_end_index == -1 and parent_end_name is None:
            raise ValueError("must supply either parent_end_index or parent_end_name")

        if parent_end_index != -1:
            return parent_end_index

        p = self._dg.get_node(parent_index)
        return p.get_end_index(parent_end_name)

    def _is_motif_in_graph_already(self, m):
        for ni in self._dg:
            n_motif = self._dg.get_node(ni)
            if m.uuid == n_motif.uuid:
                raise ValueError(
                    "cannot add motif: " + m.name + " to graph as its uuid is " +
                    "already present in the graph")
        return 0

    def add_motif(self, m, parent_index=-1, parent_end_index=0,
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
                "cannot add motif: "+ m.name + "it has steric clashes at "
                "position: " + str(self._dg.index)
            )
        # add motif to graph
        return self._dg.add_node(m_aligned, m.num_ends(), parent_index,
                                 parent_end_index, 0, index)

    def _align(self, parent_index, parent_end_index, m):
        raise NotImplementedError

    def _steric_clash(self, m):
        raise NotImplementedError

    @property
    def sterics(self):
        return self._sterics

    @sterics.setter
    def sterics(self, new_value):
        if new_value != 0 and new_value != 1:
            raise ValueError("sterics can only be set to 0 or 1")
        self._sterics = new_value

class MotifDirectedGraph(MotifTypeDirectedGraph):
    def __init__(self):
        super(self.__class__, self).__init__()
        self._motif_class_type = motif.Motif

    def nodes_to_pdbs(self, name="node"):
        for ni in self._dg:
          self._dg.get_node(ni).to_pdb(name+"."+str(ni)+".pdb")


    def _align(self, parent_index, parent_end_index, m):
        m_copy = motif.Motif.copy(m)
        if parent_index == -1:
            return m_copy

        parent = self._dg.get_parent(parent_index)
        m_copy = motif.get_aligned_motif(parent.get_end(parent_end_index), m.get_end(0), m)
        return m_copy

    def _steric_clash(self, m):
        if not self._sterics:
            return 0

        for ni in self._dg:
            if motif.clash_between_motifs(self._dg.get_node(ni), m):
                return 1
        return 0


