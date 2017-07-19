from . import exceptions

class DirectedGraph(object):
    __slots__ = [
        "_nodes",
        "_edges",
        "_level",
        "_levels",
        "_index",
        "_iter_list",
        "_parent",
        "_rebuild_list",
        "_last"]

    def __init__(self):
        self._nodes = {}
        self._edges = {}
        self._parent = {}
        self._level = 0
        self._index = 0
        self._rebuild_list = 0
        self._iter_list = []
        self._last = -1

    def __len__(self):
        return len(self._nodes)

    def __iter__(self):
        if self._rebuild_list:
            self._rebuild_iter_list()
        return self._iter_list.__iter__()

    def iter_from_root(self, ni):
        self._rebuild_list = 1
        self._rebuild_iter_list(ni)
        self._rebuild_list = 1
        return self._iter_list.__iter__()

    def get_roots(self):
        n_indexes = []
        for i, n in self._nodes.iteritems():
            if i not in self._parent:
                n_indexes.append(i)
        return n_indexes

    def get_node(self, index):
        """
        :param index: the node index that you want
        :type index: int
        :return: GraphNode object
        :raises: exceptions.GraphIndexException

        :examples:

        .. code-block:: python

            >>> g = Graph()
            >>> g.add_data(10)
            #get node of index '0' which is the first one
            >>> print g.get_node(0).data
            10
        """
        if index in self._nodes:
            return self._nodes[index]

        raise exceptions.GraphIndexException(
            index,
           "node %d was requested but does not exist" % (index))

    def get_edges(self, index):
        if index not in self._edges:
            raise ValueError(str(index) + " index does not exist in directed_graph")
        else:
            return self._edges[index]

    def get_num_edges(self, index):
        return len(self.get_edges(index))

    def get_all_edges(self):
        all_edges = []
        for ni, edges in self._edges.iteritems():
            for ei, e in enumerate(edges):
                nj, ej = e
                if nj is None:
                    continue
                full_edge = [ni, nj, ei, ej]
                if ni > e:
                    full_edge = [nj, ni, ej, ei]
                if full_edge not in all_edges:
                    all_edges.append(full_edge)
        return all_edges

    def get_parent_index(self, index):
        if index in self._parent:
            return self._parent[index]
        else:
            return None

    def get_parent(self, index):
        pi = self.get_parent_index(index)
        if pi is None:
            raise ValueError("cannot parent, node: " +
                             str(index) + " does not have a parent")
        return self.get_node(pi)

    def has_parent(self, index):
        if self.get_parent_index(index) is not None:
            return True
        else:
            return False

    def get_parent_edge_index(self, index):
        pi = self.get_parent_index(index)
        if pi is None:
            raise ValueError("cannot get parent edge index, node: " +
                             str(index) + " does not have a parent" )
        edges = self.get_edges(pi)
        for i, n in enumerate(edges):
            if n[0] == index:
                return i

    def edge_index_empty(self, ni, ei):
        edges = self.get_edges(ni)
        if ei >= len(edges):
            raise ValueError(
                "node: " + str(ni) + " has only " + str(len(edges)) + " but requested " + str(ei))

        if edges[ei][0] is None:
            return True
        else:
            return False

    def edge_exists(self, ni, nj, ei, ej):
        if self._edges[ni][ei] == [nj, ej]:
            return True
        else:
            return False

    def get_edge_index_to_parent(self, index):
        pi = self.get_parent_index(index)
        if pi is None:
            raise ValueError("cannot get edge to parent, node: " +
                             str(index) + " does not have a parent")
        edges = self.get_edges(index)
        for i, n in enumerate(edges):
            if pi == n[0]:
                return i

    def are_nodes_connected(self, ni, nj):
        edges = self.get_edges(ni)

        for connected_n, connected_end in edges:
            if connected_n == nj:
                return True

        return False

    def get_connected_node(self, ni, ei):
        return self._edges[ni][ei]

    def add_edge(self, ni, nj, ei, ej):
        if not self.edge_index_empty(ni, ei):
            raise ValueError(
                "attempting to add edge to node: " + str(ni) +
                " and edge_index: " + str(ei) + " is already filled by node " +
                str(self.get_connected_node(ni, ei)))

        if not self.edge_index_empty(nj, ej):
            raise ValueError(
                "attempting to add edge to node: " + str(nj) +
                " and edge_index: " + str(ej) + " is already filled by node " +
                str(self.get_connected_node(nj, ej)))

        self._edges[ni][ei] = [nj, ej]
        self._edges[nj][ej] = [ni, ei]
        self._rebuild_list = 1

    def add_node(self, n, allowed_edges, parent_index=-1, parent_edge_index=1,
                 n_edge_index=0, index=-1):

        # what is the index of this node?
        i = index
        if i != -1:
            self._check_if_unused_index(i)
            self._index = i
        else:
            i = self._index

        self._nodes[i] = n
        self._edges[i] = [[None, None] for j in range(allowed_edges)]

        self._index = self._get_next_available_index()

        if parent_index != -1:
            self.add_edge(i, parent_index, n_edge_index, parent_edge_index)
            self._parent[i] = parent_index

        self._rebuild_list = 1
        self._last = i

        return i

    def add_graph(self, dg, node_index, node_edge_index,
                  parent_index=-1, parent_edge_index=1):

        updated_indexes = {}
        first = 1
        for ni in dg.iter_from_root(node_index):
            n = dg.get_node(ni)
            if dg.has_parent(ni):
                new_parent_index = updated_indexes[ dg.get_parent_index(ni) ]
                new_index = self.add_node(
                    n, dg.get_num_edges(ni), new_parent_index,
                    dg.get_parent_edge_index(ni), dg.get_edge_index_to_parent(ni))
            else:
                if first:
                    first = 0
                    new_index = self.add_node(
                        n, dg.get_num_edges(ni), parent_index, parent_edge_index,
                        node_edge_index)
                else:
                    new_index = self.add_node(n, dg.get_num_edges(ni))
            updated_indexes[ni] = new_index

        # catch extra connections not defined from parent to child
        for full_edge in dg.get_all_edges():
            ni = updated_indexes [ full_edge[0] ]
            nj = updated_indexes [ full_edge[1] ]
            if not self.are_nodes_connected(ni, nj):
                self.add_edge(ni, nj, full_edge[2], full_edge[3])

    def remove_node(self, ni):
        if ni not in self._nodes:
            raise ValueError("cannot remove node: " + str(ni) + " it does not exist")

        edges = self.get_edges(ni)
        for i, n in enumerate(edges):
            nj, ej = n
            if nj is None:
                continue
            # no longer has a parent
            if nj in self._parent:
                if self._parent[nj] == ni:
                    del self._parent[nj]

            self._edges[nj][ej] = [None, None]


        del self._edges[ni]
        if ni in self._parent:
            del self._parent[ni]
        del self._nodes[ni]

        self._rebuild_list = 1

    def remove_edge(self, ni, nj, ei, ej):
        if self._edges[ni][ei] != [nj, ej]:
            raise ValueError(
                "edge: %d %d %d %d" % (ni, nj, ei, ej) + "does not exist" )

        if ni == self.get_parent_index(nj):
            del self._parent[nj]

        if nj == self.get_parent_index(ni):
            del self._parent[ni]

        self._edges[nj][ej] = [None, None]
        self._edges[ni][ei] = [None, None]
        self._rebuild_list = 1

    def _rebuild_iter_list(self, ni=None):
        self._iter_list = []
        roots = self.get_roots()
        open = []
        seen = []
        for root in roots:
            open.append(root)
            seen.append(root)

            while len(open) > 0:
                ni = open.pop(0)

                for i in self._nodes.iterkeys():
                    if i not in self._parent:
                        continue
                    if self._parent[i] == ni and i not in seen:
                        open.append(i)
                        seen.append(i)

                self._iter_list.append(ni)
        self._rebuild_list = 0

    def _check_if_unused_index(self, index):
        if index in self._nodes:
            raise exceptions.GraphIndexException(
                index,
                "specified index: %d cannot be used as it already" % (index) + \
                " assigned to a node")
        else:
            return 1

    def _get_next_available_index(self):
        next_index = self._index + 1
        while 1:
            if next_index not in self._nodes:
                return next_index
            else:
                next_index += 1

    # getters
    @property
    def index(self):
        return self._index

    @property
    def last(self):
        return self._last




