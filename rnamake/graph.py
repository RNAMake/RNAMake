import Queue

def transverse_graph(graph, i):
    graph.__iter__()
    graph.current_node = graph.get_node(i)
    while 1:
        try:
            next = graph.next()
        except:
            raise StopIteration
        yield next


class Graph(object):
    """
    A general graph storage base class. Do not call directly!

    Attributes
    ----------
    `nodes` : list of GraphNode objects
        all the nodes in the tree, used for fast access by index
    `connections` : list of GraphConnection objects
        all the connections between two different nodes in the graph
    `level` : int
        the current graph level, used for quickly deleting sections of the graph
    `index`: int
        the current node index, always the length of the number of nodes in teh graph
    `last_node` : GraphNode object
        the last node added to the graph
    `current_node` GraphNode object
        the current node during iteration

    """

    def __init__(self):
        self.nodes, self.connections, self.level, self.index = [], [], 0, 0
        self.last_node = None
        #iterator stuff
        self.current_node = None

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        self.current_node = self.nodes[0]
        return self

    def next(self):
        if self.current_node is None:
            raise StopIteration

        node = self.current_node

        try:
            self.current_node = self.get_node(self.current_node.index+1)
        except:
            self.current_node = None

        return node

    def get_node(self, index):
        """
        :param index: the node index that you want
        :type index: int
        :return: GraphNode object

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

        raise ValueError("cannot find node with index")


class GraphDynamic(Graph):
    """
    a Graph with dynamic connections between nodes. i.e. each node does NOT have a predefined
    number of connections. Each node starts with 0 connections and are added over time,
    there is no max to the number of connections each node can have

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

        :param data:
        :param parent_index:
        :return:
        """
        parent = None
        if parent_index != -1:
            parent = self.get_node(parent_index)

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
        n1 = self.get_node(i)
        n2 = self.get_node(j)
        c = GraphConnection(n1, n2, 0, 0)
        n1.add_connection(c)
        if n1 != n2:
            n2.add_connection(c)
        self.connections.append(c)


class GraphStatic(Graph):
    def __init__(self):
        super(GraphStatic, self).__init__()

    def add_data(self, data, parent_index=-1, parent_pos=-1, child_pos=-1, n_children=1):
        parent = self.last_node
        if parent_index != -1:
            parent = self.get_node(parent_index)
        n = GraphNodeStatic(data, self.index, self.level, n_children)

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
        return self.index-1

    def connect(self, i, j, i_pos, j_pos):
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
        self.last_node = self.nodes[-1]
        #self.index -= 1

    def check_pos_is_value(self, n, pos):
        if pos == -1:
            avail_pos = n.available_children_pos()
            if len(avail_pos) == 0:
                raise ValueError("cannot add connection to node has no available ends")
            return avail_pos[0]

        else:
            if n.available_pos(pos) == 0:
                raise ValueError("graph pos is not available")
            return pos

    def get_availiable_pos(self, n, pos):
        if pos == -1:
            avail_pos = n.available_children_pos()
            return avail_pos
        else:
            if n.available_pos(pos) == 0:
                raise ValueError("graph pos is not available: " + str(pos) + " " + n.data.name)
            return [pos]



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
