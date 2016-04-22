import exceptions
import priority_queue


# TODO this is actually not a great system should probably be a member of Graph
def transverse_graph(graph, i, directed=-1):
    graph.__iter__()
    graph.current_node = graph.get_node(i)
    graph.directed = 0
    while 1:
        try:
            next = graph.next()
        except:
            raise StopIteration
        yield next


# TODO using a priority queue is not enough to assume connectivity will be
# followed correctly. Need to inforce that connectivity will go from
# motif that is aligned to one that is not and not build connections between
# ends 1 of two different motifs such at the end of a build path possibly need
# seperate class that is a DirectedGraph instead of just Graph
# see /rnamake/unittests/resources/motif_graph/test_added.mg to see this issue
class Graph(object):
    """
    General implementation of a undirected graph. Do not call directly!

    :attributes:
    `nodes` : list of GraphNode objects
        all the nodes in the tree, used for fast access by index
    `connections` : list of GraphConnection objects
        all the connections between two different nodes in the graph
    `level` : int
        the current graph level, used for quickly deleting sections of the graph
    `index`: int
        the current node index, always the length of the number of nodes in teh graph
    `last_node`: GraphNode object
        the last node added to the graph
    `current_node`: GraphNode object
        the current node during iteration
    `queue`: PriorityQueue object
        tracks which nodes to visit turning iteration
    `seen`: List of GraphNode Objects
        tracks which nodes have alaready been visited during iteration
    """

    def __init__(self):
        self.nodes, self.connections, self.level, self.index = [], [], 0, 0
        self.last_node = None
        # iterator stuff
        self.current_node = None
        self.queue = priority_queue.PriorityQueue()
        self.seen = []
        self.directed = -1

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        if len(self.nodes) != 0:
            for n in self.nodes:
                avail_pos = n.available_children_pos()
                if len(avail_pos) != 0:
                    self.current_node = n
                    break
            if self.current_node is None:
                self.current_node = self.nodes[0]
        else:
            self.current_node = None
        self.seen = [self.current_node]
        self.queue = priority_queue.PriorityQueue()
        self.directed = -1
        return self

    def next(self):
        if self.current_node is None:
            raise StopIteration

        node = self.current_node

        for i, c in enumerate(node.connections):
            if c is None:
                continue

            partner = c.partner(node.index)
            if c.end_index(partner.index) != 0 and self.directed == 0:
                continue
            if partner not in self.seen:
                self.queue.put(partner, partner.index)
                self.seen.append(partner)

        if self.queue.empty() and len(self.seen) == len(self.nodes):
            self.current_node = None
        elif self.queue.empty():
            self.current_node = None

            for n in self.nodes:
                if n not in self.seen:
                    self.seen.append(n)
                    self.current_node = n
                    break

        else:
            self.current_node = self.queue.get()


        return node

    def get_node(self, index):
        """
        :param index: the node index that you want
        :type index: int
        :return: GraphNode object

        :examples:

        .. code-block:: python

            >>>g = Graph()
            >>>g.add_data(10)
            #get node of index '0' which is the first one
            >>>print g.get_node(0).data
            10
        """
        for n in self.nodes:
            if n.index == index:
                return n

        raise exceptions.GraphException("cannot find node with index: " + str(index))

    def oldest_node(self):
        """
        returns the node with the lowest index. Only used a few examples may
        deprecated in next verision.

        :returns: GraphNode object

        :examples:

        .. code-block:: python

            >>>g = Graph()
            >>>g.add_data(10)
            >>>g.add_data(5)
            >>>g.add_data(15)
            >>> g.oldest_node().data
            10
        """
        node = self.last_node
        for n in self.nodes:
            if n.index < node.index:
                node = n
        return node

    def increase_level(self):
        """
        Increases the level of nodes to be added. default level is 0. This is
        useful when removing or adding a set of nodes. Think of level as a
        grouping mechanism
        """
        self.level += 1

    def decrease_level(self):
        """
        Decreases the level of nodes to be added. default level is 0. This is
        useful when removing or adding a set of nodes. Think of level as a
        grouping mechanism
        """

        if self.level == 0:
            raise exceptions.GraphException("cannot decrease level anymore is "
                                            "already 0")

        self.level -= 1


# TODO this could probably just be Graph doesnt need to be its own class
class GraphDynamic(Graph):
    """
    a Graph with dynamic connections between nodes. i.e. each node does NOT
    have a predefined number of connections. Each node starts with 0
    connections and are added over time, there is no max to the number of
    connections each node can have

    Examples

    .. code-block:: python

        >>>g = GraphDynamic()
        >>>g.add_data(0)
        >>>g.add_data(1)
        >>>g.add_data(2, parent_index=0)
        >>>g.add_data(3, parent_index=0)

        >>>len(g.get_node(0).connections)
        3
    """

    def __init__(self):
        super(GraphDynamic, self).__init__()

    def add_data(self, data, parent_index=-1):
        """
        add a new peice of data to the graph

        :param data: Item to be added to graph
        :param parent_index: index of node that this element should be connected too

        :type data: can be anything
        :type parent_index: int

        :return: index of added node
        :rtype: int
        """
        parent = None
        if parent_index != -1:
            try:
                parent = self.get_node(parent_index)
            # cannot find node of that index
            except exceptions.GraphException:
                raise exceptions.GraphException("cannot add data to parent_" +
                                                "index: " + str(parent_index) +
                                                " that nodes does not exist")

        n = GraphNodeDynamic(data, self.index, self.level)

        if parent is not None:
            c = GraphConnection(parent, n, 0, 0)
            parent.add_connection(c)
            n.add_connection(c)
            self.connections.append(c)

        self.nodes.append(n)
        self.index += 1
        self.last_node = n
        return self.index - 1

    def connect(self, i, j):
        """
        Add a connection between two nodes. Since this is a Dynamic Graph this
        cannot fail. i and j can be the same

        :param i: index of node 1 to connect
        :param j: index of node 2 to connect

        :type i: int
        :type j: int
        :return: None
        """
        n1 = self.get_node(i)
        n2 = self.get_node(j)
        c = GraphConnection(n1, n2, 0, 0)
        n1.add_connection(c)
        if n1 != n2:
            n2.add_connection(c)
        self.connections.append(c)


class GraphStatic(Graph):
    """
    a Graph with static connections between nodes. i.e. each node HAS a
    predefined number of connections. This is useful for constructing graphs
    of Motifs and Secondary Structure elements which have predefined ends.

    :examples:

    .. code-block:: python

        >>> g = graph.GraphStatic()
        # legend
        # N0 = Node 0
        # N0C0 = Connection 0 of Node 0
        # add first node with 2 possible connection
        >>> g.add_data(0, n_children=2)
        #          N0
        #    N0C0 /  \  N0C1
        #        /    \\
        #     None     None
        # add new node, connected to node 0. This new node: node 1 is connected
        # to node 0 using the connnection in the 0th position on both node 0 and
        # node 1. Do not need to specifiy parent_index=0 since adds to last node
        # without specifying
        # g.add_data(1, parent_index=0, parent_pos=0, child_pos=0, n_children=2)
        # means the same thing
        >>> g.add_data(1, parent_pos=0, child_pos=0, n_children=2)
        #          N0
        #    N0C0 /  \  N0C1
        #    N1C0/    \\
        #       N1   None
        #   N1C1|
        #       |
        #      None
        # add new node, connected to node 0.
        >>> g.add_data(2, parent_index=0, parent_pos=1, child_pos=0, n_children=2)
        #          N0
        #    N0C0 /  \  N0C1
        #    N1C0/    \ N2C0
        #       N1    N2
        #   N1C1|      |N2C1
        #       |      |
        #      None   None
        >>> g.connect(1, 2, 1, 1)
        #          N0
        #    N0C0 /  \  N0C1
        #    N1C0/    \ N2C0
        #       N1----N2
        #      N1C1  N2C1
    """

    def __init__(self):
        super(GraphStatic, self).__init__()

    def add_data(self, data, parent_index=-1, parent_pos=-1, child_pos=-1, n_children=1,
                 orphan=0, index=-1):
        """
        Adds a new node to the graph given specfied data.

        :param data: element to add to graph
        :param parent_index: index of parent to connect to, optional
        :param parent_pos:  connection position to connect to parent, optional
        :param child_pos:  connection position on new node to be created, optional
        :param n_children: number of connections new node will have, optional
        :param orphan: whether node will not be connected to any other node.
            default new nodes are connected to the last node added. if orphan=1
            this will not happen
        :param index: overrides internal index and set node to specific index
            value, this is rarely used other then copying an entire graph

        :type parent_index: int
        :type parent_pos: int
        :type child_pos: int
        :type n_children: int
        :type orphan: int
        :type index: int

        :return: index of new node
        :rtype: int

        :examples:

        .. code-block:: python

            >>> g = graph.GraphStatic()
            # legend
            # N0 = Node 0
            # N0C0 = Connection 0 of Node 0
            # add first node with 2 possible connection
            >>> g.add_data(0, n_children=2)
            #          N0
            #    N0C0 /  \  N0C1
            #        /    \\
            #     None     None
            # add new node, connected to node 0. This new node: node 1 is connected
            # to node 0 using the connnection in the 0th position on both node 0 and
            # node 1. Do not need to specifiy parent_index=0 since adds to last node
            # without specifying
            # g.add_data(1, parent_index=0, parent_pos=0, child_pos=0, n_children=2)
            # means the same thing
            >>> g.add_data(1, parent_pos=0, child_pos=0, n_children=2)
            #          N0
            #    N0C0 /  \  N0C1
            #    N1C0/    \\
            #       N1   None
            #   N1C1|
            #       |
            #      None
        """

        given_index = self.index
        if index != -1:
            given_index = index
        parent = self.last_node
        if parent_index != -1:
            parent = self.get_node(parent_index)
        if orphan:
            parent = None

        n = GraphNodeStatic(data, given_index, self.level, n_children)

        if parent is not None:
            parent_pos = self.check_pos_is_value(parent, parent_pos)
            child_pos  = self.check_pos_is_value(n, child_pos)
            c = GraphConnection(parent, n, parent_pos, child_pos)
            parent.add_connection(c, parent_pos)
            n.add_connection(c, child_pos)
            self.connections.append(c)

        self.nodes.append(n)
        self.index += 1
        self.last_node = n
        return given_index

    def connect(self, i, j, i_pos, j_pos):
        """
        Connects two nodes together

        :param i:
        :param j:
        :param i_pos:
        :param j_pos:
        :return:
        """
        n1 = self.get_node(i)
        n2 = self.get_node(j)
        i_pos = self.check_pos_is_value(n1, i_pos)
        j_pos = self.check_pos_is_value(n2, j_pos)
        c = GraphConnection(n1, n2, i_pos, j_pos)
        n1.add_connection(c, i_pos)
        n2.add_connection(c, j_pos)

        self.connections.append(c)

    def remove_node(self, pos):
        n = self.get_node(pos)
        for c in n.connections:
            if c is not None:
                c.disconnect()
                self.connections.remove(c)
        self.nodes.remove(n)
        self.last_node = None
        if len(self.nodes) > 0:
            self.last_node = self.nodes[-1]
        #self.index -= 1

    def remove_node_level(self, level=None):
        if level is None:
            level = self.level

        nodes = self.nodes[::-1]
        while 1:
            removed = 0
            for n in nodes:
                if n.level == level:
                    self.remove_node(n.index)
                    removed = 1
                    break
            nodes = self.nodes[::-1]
            if not removed:
                break

    def check_pos_is_value(self, n, pos, error=1):
        if pos == -1:
            avail_pos = n.available_children_pos()
            if len(avail_pos) == 0:
                raise ValueError("cannot add connection to node has no available ends")
            return avail_pos[0]

        else:
            if n.available_pos(pos) == 0 and error:
                raise ValueError("graph pos is not available: " + str(pos) +
                                 " on node: " + str(n.index) )
            elif n.available_pos(pos) == 0:
                return None
            return pos

    def get_availiable_pos(self, n, pos):
        if pos == -1:
            avail_pos = n.available_children_pos()
            return avail_pos
        else:
            if n.available_pos(pos) == 0:
                raise ValueError("graph pos is not available: " + str(pos) + " " + n.data.name)
            return [pos]

    def copy(self):
        gs = GraphStatic()
        new_nodes = []
        for n in self.nodes:
            new_nodes.append(n.copy())
        gs.nodes = new_nodes
        for c in self.connections:
            i, j = c.node_1.index, c.node_2.index
            ni, nj = c.end_index_1, c.end_index_2
            gs.connect(i, j, ni, nj)
        if self.last_node is not None:
            gs.last_node = gs.nodes[self.last_node.index]
        gs.level = self.level
        gs.index = self.index
        return gs


class GraphNode(object):
    def __init__(self, data, index, level, n_connections=0):
        self.data, self.index, self.level = data, index, level
        self.connections = [None for x in range(n_connections)]

    def available_children_pos(self):
        pos = []
        for i, c in enumerate(self.connections):
            if c is None:
                pos.append(i)
        return pos

    def available_pos(self, pos):
        if len(self.connections) <= pos:
            return 0
        if self.connections[pos] is not None:
            return 0
        return 1

    def parent(self):
        for c in self.connections:
            if c is None:
                continue
            n = c.partner(self.index)
            if n.index < self.index:
                return n

        return None

    def parent_index(self):
        parent = self.parent()
        if parent is None:
            return -1
        for c in self.connections:
            if c.partner(self.index) == parent.index:
                return c.end_index(parent.index)

        return -1

    def connected(self, n):
        for c in self.connections:
            if c is None:
                continue
            if c.partner(self.index) == n:
                return c
        return None


class GraphNodeDynamic(GraphNode):
    def __init__(self, data, index, level):
        super(GraphNodeDynamic, self).__init__(data, index, level)

    def add_connection(self, c):
        self.connections.append(c)

    def remove_connection(self, c):
        if c not in self.connections:
            raise ValueError("cannot remove connection not in node")

        self.connections.remove(c)


class GraphNodeStatic(GraphNode):
    def __init__(self, data, index, level, n_children):
        super(GraphNodeStatic, self).__init__(data, index, level, n_children)

    def add_connection(self, c, pos):
        if pos >= len(self.connections):
            raise ValueError("cannot add child at position")

        if self.connections[pos] is not None:
            raise ValueError("attempted to add child in a position that is already full")

        self.connections[pos] = c

    def remove_connection(self, c):
        if c not in self.connections:
            raise ValueError("cannot remove connection not in node")

        i = self.connections.index(c)
        self.connections[i] = None

    def copy(self):
        new_data = None
        try:
            new_data = self.data.copy()
        except:
            new_data = self.data

        c = GraphNodeStatic(new_data, self.index, self.level, len(self.connections))
        return c


#TODO rename connection position from end. Not general enough
class GraphConnection(object):
    def __init__(self, node_1, node_2, end_index_1, end_index_2):
        self.node_1, self.node_2 = node_1, node_2
        self.end_index_1, self.end_index_2 = end_index_1, end_index_2

    def partner(self, n_index):
        if   n_index == self.node_1.index:
            return self.node_2
        elif n_index == self.node_2.index:
            return self.node_1
        else:
            raise ValueError("cannot call partner with node not in connection, index: " + \
                              n_index + " if this is not a number you messed up")

    def end_index(self, n_index):
        if   n_index == self.node_1.index:
            return self.end_index_1
        elif n_index == self.node_2.index:
            return self.end_index_2
        else:
            raise ValueError("cannot call end_index with node not in connection")

    def disconnect(self):
        self.node_1.remove_connection(self)
        self.node_2.remove_connection(self)
