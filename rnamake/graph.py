import Queue

class Graph(object):
    def __init__(self):
        self.nodes, self.connections, self.level, self.index = [], [], 0, 0
        self.last_node = None
        #iterator stuff
        self.current_node = None
        self.queue = Queue.Queue()
        self.seen = {}

    def get_node(self, index):
        for n in self.nodes:
            if n.index == index:
                return n

        raise ValueError("cannot find node with index")

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        self.current_node = self.get_node(0)
        while not self.queue.empty():
            self.queue.get()
        self.seen = {}
        self.seen[self.current_node] = 1
        return self

    def next(self):
        if self.current_node is None:
            raise StopIteration

        for c in self.current_node.connections:
            if c is None:
                continue
            partner = c.partner(self.current_node.index)
            if partner in self.seen:
                continue

            self.seen[partner] = 1
            self.queue.put(partner)

        node = self.current_node
        if self.queue.empty():
            self.current_node = None
        else:
            self.current_node = self.queue.get()

        return node


class GraphDynamic(Graph):
    def __init__(self):
        super(GraphDynamic, self).__init__()

    def add_data(self, data, parent_index=-1):
        parent = self.last_node
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
            parent_pos = self._check_pos_is_value(parent, parent_pos)
            child_pos  = self._check_pos_is_value(n, child_pos)
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
        i_pos = self._check_pos_is_value(n1, i_pos)
        j_pos = self._check_pos_is_value(n2, j_pos)
        c = GraphConnection(n1, n2, i_pos, j_pos)
        n1.add_connection(c, i_pos)
        n2.add_connection(c, j_pos)
        self.connections.append(c)

    def remove_node(self, pos):
        n = self.get_node(pos)
        for c in n.connections:
            c.disconnect()
            self.connections.remove(c)
        self.nodes.remove(n)


    def _check_pos_is_value(self, n, pos):
        if pos == -1:
            avail_pos = n.available_children_pos()
            if len(avail_pos) == 0:
                raise ValueError("cannot add connection to node has no available ends")
            return avail_pos[0]

        else:
            if n.available_pos(pos) == 0:
                raise ValueError("graph pos is not available")
            return pos


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

    def n_connections(self):
        count = 0
        for c in self.connections:
            if c is not None:
                count += 1
        return count


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
            raise ValueError("cannot call partner with not not in connection")

    def disconnect(self):
        self.node_1.remove_connection(self)
        self.node_2.remove_connection(self)
