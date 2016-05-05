import exceptions

# TODO should be moved into Tree
def transverse_tree(tree, i):
    tree.__iter__()
    tree.current_node = tree.get_node(i)
    while 1:
        try:
            next = tree.next()
        except:
            raise StopIteration
        yield next

class Tree(object):
    """
    General implementation of a tree. Do not call directly!

    :attributes:
    `nodes` : list of TreeNode objects
        all the nodes in the tree, used for fast access by index
    `level` : int
        the current tree level, used for quickly deleting sections of the tree
    `index`: int
        the current node index, always the length of the number of nodes in the tree
    `last_node`: TreeNode object
        the last node added to the tree
    `current_node`: TreeNode object
        the current node during iteration
    """

    def __init__(self):
        self.nodes = []
        self.last_node = None
        self.index = 0
        self.level = 0

        #iterator stuff
        self.current_node = None

    def get_node(self, index):
        """
        :param index: the node index that you want
        :type index: int
        :return: TreeNode object

        :examples:

        .. code-block:: python

            >>> t = TreeDynamic()
            >>> t.add_data(10)
            #get node of index '0' which is the first one
            >>> print t.get_node(0).data
            10
        """

        for n in self.nodes:
            if n.index == index:
                return n

        raise exceptions.TreeException(
            "cannot find node with index "+ str(index))

    def remove_node(self, node=None, index=None):
        """
        removes a node with a given index

        :param node: the actual node object you want to remove
        :param index: the index of the node you want to remove

        :type node: TreeNode
        :type index: int

        :return: None
        """

        if node is None and index is not None:
            node = self.get_node(index)
        elif node is None and index is None:
            raise exceptions.TreeException(
                "need to include node or index in remove_node")

        if node.parent is not None:
            node.parent.remove_child(node)
        self.nodes.remove(node)
        node.parent = None

        if self.last_node == node:
            self.last_node = node.parent

    def remove_node_level(self, level=None):
        """
        remove all nodes with a given level, this is useful if you are unsure
        how many nodes have been added from a given point. if level is not
        specified level will be set to current level

        :param level: the node minimum node level you want to remove
        :type level: int
        :return: None

        """

        if level is None:
            level = self.level

        nodes = self.nodes[::-1]
        while 1:
            removed = 0
            for n in nodes:
                if n.level == level:
                    self.remove_node(n)
                    removed = 1
                    break
            nodes = self.nodes[::-1]
            if not removed:
                break

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        self.current_node = self.get_node(0)
        return self

    # TODO this is really bad need to just call connections.__iter__()
    def next(self):
        if self.current_node is None:
            raise StopIteration

        node = self.current_node

        if len(self.nodes)-1 == node.index:
            self.current_node = None
        else:
            self.current_node = self.nodes[node.index + 1]

        return node

    def next_level(self):
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

        self.level -= 1


class TreeDynamic(Tree):
    """
    a Tree that has a dynamic number of children. i.e. each node does NOT
    have a predefined number of children. Each node starts with 0
    children and are added over time, there is no max to the number of
    children each node can have

    :attributes:
    `nodes` : list of TreeNode objects
        all the nodes in the tree, used for fast access by index
    `level` : int
        the current tree level, used for quickly deleting sections of the tree
    `index`: int
        the current node index, always the length of the number of nodes in the tree
    `last_node`: TreeNode object
        the last node added to the tree
    `current_node`: TreeNode object
        the current node during iteration

    .. code-block:: python

        >>> t = TreeDynamic()
        >>> t.add_data(0)
        >>> t.add_data(1)
        >>> t.add_data(2, parent_index=0)
        >>> t.add_data(3, parent_index=0)

        >>> len(t.get_node(0).children)
        3

    """

    def __init__(self):
        super(TreeDynamic, self).__init__()

    def add_data(self, data, parent_index=-1):
        """
        add a new element of data to the tree

        :param data: Item to be added to tree
        :param parent_index: index of node that this element should be connected too

        :type data: can be anything
        :type parent_index: int

        :return: index of added node
        :rtype: int
        """

        n = TreeNodeDynamic(data, self.index, self.level)
        parent = self.last_node
        if parent_index != -1:
            try:
                parent = self.get_node(parent_index)
            except exceptions.TreeException:
                raise exceptions.TreeException(
                    "cannot add data to parent index: " + str(parent_index) +
                    " that nodes does not exist")

        if parent != None:
            parent.add_child(n)
            n.parent = parent

        self.last_node = n
        self.nodes.append(n)
        self.index += 1
        return self.index-1


class TreeStatic(Tree):
    def __init__(self):
        super(TreeStatic, self).__init__()

    def add_data(self, data, n_children=1, parent_index=-1, parent_child_index=-1):
        n = TreeNodeStatic(data, self.index, self.level, n_children)
        parent = self.last_node

        if parent_index != -1:
            parent = self.get_node(parent_index)

        if parent != None:
            if parent_child_index == -1:
                all_pos = parent.available_children_pos()
                if len(all_pos) == 0:
                    raise ValueError("cannot add child to this parent it has no open spots")
                parent_child_index = all_pos[0]
            parent.add_child(n, parent_child_index)
            n.parent = parent

        self.last_node = n
        self.nodes.append(n)
        self.index += 1
        return self.index-1

    def get_available_pos(self, n, pos):
        if pos == -1:
            avail_pos = n.available_children_pos()
            return avail_pos
        else:
            if n.available_pos(pos) == 0:
                raise ValueError("tree pos is not available")
            return [pos]

    def copy(self):
        ts = TreeStatic()
        new_nodes = []
        for n in self.nodes:
            new_nodes.append(n.copy())
        for n in self.nodes:
            parent_index = n.parent_index()
            pei = n.parent_end_index()
            if parent_index == -1:
                continue
            new_nodes[parent_index].add_child(n, pei)
            n.parent = new_nodes[parent_index]

        ts.nodes = new_nodes
        if self.last_node is not None:
            ts.last_node = ts.nodes[self.last_node.index]
        ts.level = self.level
        ts.index = self.index
        return ts



class TreeNode(object):
    def __init__(self, data, index, level, n_children=0):
        self.parent = None
        self.data, self.index, self.level = data, index, level
        self.children = [ None for x in range(n_children)]

    def available_children_pos(self):
        pos = []
        i = 0
        for i, c in enumerate(self.children):
            if c is None:
                pos.append(i)

        return pos

    def parent_index(self):
        if self.parent == None:
            return -1

        else:
            return self.parent.index

    def available_pos(self, pos):
        if len(self.children) <= pos:
            return 0
        if self.children[pos] is not None:
            return 0
        return 1

    def parent_end_index(self):
        if self.parent is None:
            return -1
        return self.parent.children.index(self)

class TreeNodeDynamic(TreeNode):
    def __init__(self, data, index, level):
        super(TreeNodeDynamic, self).__init__(data, index, level)

    def add_child(self, node):
        self.children.append(node)

    def remove_child(self, node):
        if node not in self.children:
            raise ValueError("cannot remove child from node as it is not a child")

        self.children.remove(node)


class TreeNodeStatic(TreeNode):
    def __init__(self, data, index, level, n_children):
        super(TreeNodeStatic, self).__init__(data, index, level, n_children)

    def add_child(self, node, pos):
        if pos >= len(self.children):
            raise ValueError("pos is greater then children slots")

        if self.children[pos] is not None:
            raise ValueError("cannot add child at this spot there is already a child here")

        self.children[pos] = node

    def remove_child(self, node):
        if node not in self.children:
            raise ValueError("cannot remove child from node as it is not a child")

        i = self.children.index(node)
        self.children[i] = None

    def copy(self):
        new_data = None
        try:
            new_data = self.data.copy()
        except:
            new_data = self.data
        c = TreeNodeStatic(new_data, self.index, self.level, len(self.children))
        return c


