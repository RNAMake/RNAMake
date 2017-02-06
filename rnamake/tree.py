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

        # iterator stuff
        self.current_node = None

    def get_node(self, index):
        """
        :param index: the node index that you want
        :type index: int
        :return: TreeNode object
        :raises: exceptions.TreeIndexException

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

        raise exceptions.TreeIndexException(
            index,
            "node %d was requested but does not exist" % (index))

    # TODO this is not behaving properly, should be removing all children too
    # tree does not work with there are gaps in it like a graph
    def remove_node(self, node=None, index=None):
        """
        removes a node with a given index

        :param node: the actual node object you want to remove
        :param index: the index of the node you want to remove

        :type node: TreeNode
        :type index: int

        :return: None
        :raises: exceptions.TreeIndexException, exceptions.TreeException
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
        for c in node.children:
            if c is not None:
                c.parent = None

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

    # TODO this is really bad need to just call nodes.__iter__()
    def next(self):
        if self.current_node is None:
            raise StopIteration

        node = self.current_node

        index = self.nodes.index(node)

        if len(self.nodes)-1 == index:
            self.current_node = None
        else:
            self.current_node = self.nodes[index + 1]

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

    :examples:

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

    # TODO parent adding should be done inside adding a child, should not be
    # a seperate step
    def add_data(self, data, parent_index=-1):
        """
        add a new element of data to the tree

        :param data: Item to be added to tree
        :param parent_index: index of node that this element should be connected too

        :type data: can be anything
        :type parent_index: int

        :return: index of added node
        :rtype: int
        :raises: exceptions.TreeIndexException
        """

        n = TreeNodeDynamic(data, self.index, self.level)
        parent = self.last_node
        if parent_index != -1:
            parent = self.get_node(parent_index)

        if parent is not None:
            parent.add_child(n)
            n.parent = parent

        self.last_node = n
        self.nodes.append(n)
        self.index += 1
        return self.index-1


class TreeStatic(Tree):
    """
    A Tree that has a static number of children. i.e. each node has
    a predefined number of children. Each node starts with 0
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

    """

    def __init__(self):
        super(TreeStatic, self).__init__()

    # TODO decide what the correct order of where n_children should be its
    # different from this same function in graph.GraphStatic implementation
    def add_data(self, data, n_children=1, parent_index=-1, parent_child_index=-1):
        """
        Adds a new node to the tree given specfied data.

        :param data: element to add to tree
        :param n_children: number of children the node can have in max
        :param parent_index: index of parent to connect to, optional
        :param parent_child_index:  child position to connect to parent, optional

        :type data: anything
        :type n_children: int
        :type parent_index: int
        :type parent_child_index: int

        :return: index of new node
        :rtype: int
        :raises: exceptions.TreeIndexException
        """

        n = TreeNodeStatic(data, self.index, self.level, n_children)
        parent = self.last_node

        if parent_index != -1:
            parent = self.get_node(parent_index)

        if parent is not None:
            if parent_child_index == -1:
                all_pos = parent.available_children_pos()
                if len(all_pos) == 0:
                    raise exceptions.TreeEndIndexException(
                        parent,
                        parent_child_index,
                        """cannot add new node to parent %d since it has not
                           available ends to add children""" % (parent.index))

                parent_child_index = all_pos[0]

            parent.add_child(n, parent_child_index)
            n.parent = parent

        self.last_node = n
        self.nodes.append(n)
        self.index += 1
        return self.index-1

    # TODO get rid of this function is useless and bad
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
    """
    An abstract class to define the behavior of a tree node or an element of
    data that has connectivity relationships to other elements of data of the
    same type. This class is never called directly.

    :param data: The element that you wish to store in the tree
    :param index: The way for you to find this element of data again through
        :func:`Tree.get_node`
    :param level: A way of grouping nodes together, increasing the level and
        then removing all nodes of a given level is a quick way to batch
        remove nodes
    :param n_children: The number of connections this node can have, in
        TreeDynamic this does nothing.

    :type data: anything
    :type index: int
    :type level: int
    :type n_children: int

    :attributes:

    `data` : anything
        The element of data being stored in tree
    `index`: int
        The way for you to find this element of data again through
        :func:`Tree.get_node`
    `level`: int
        A way of grouping nodes together, increasing the level and
        then removing all nodes of a given level is a quick way to batch
        remove nodes
    `children`: list of TreeNodes
        Nodes whose parent is this node
    `parent`: TreeNode
        the parent of the current node

    """

    __slots__ = ["data", "index", "level", "children", "parent"]

    def __init__(self, data, index, level, n_children=0):
        self.parent = None
        self.data, self.index, self.level = data, index, level
        self.children = [ None for x in range(n_children)]

    def available_children_pos(self):
        """
        gets all of the children positions that are not filled yet for this
        node

        :return: a list of end positions not yet filled
        :rtype: list of ints
        """
        pos = []
        i = 0
        for i, c in enumerate(self.children):
            if c is None:
                pos.append(i)

        return pos

    def parent_index(self):
        """
        gets the node index of parent if parent is None returns -1

        :return: node index of parent
        :rtype: int
        """

        if self.parent is None:
            return -1
        else:
            return self.parent.index

    def available_pos(self, pos):
        """
        checks to see if a child position is filled or not. returns 0 is filled
        or not available and 1 if its avialable to be used.

        :param pos: the position in children list you want to check
        :return: the avialability of the position
        :rtype: int
        """

        if len(self.children) <= pos:
            return 0
        if self.children[pos] is not None:
            return 0
        return 1

    def parent_end_index(self):
        """
        gets the position of this node in the children list of the parent.
        i.e. what number child is this not in reference to the parent. will
        return -1 if this node has no parent.

        :return: the child position this node is in reference to its parent
        :rtype: int
        """

        if self.parent is None:
            return -1
        return self.parent.children.index(self)


class TreeNodeDynamic(TreeNode):
    """
    A TreeNode container that is specific for TreeDynamic trees. That is
    it can be connected to an unlimited number of other nodes. There is no
    reason to instantiate this class outside TreeDynamic.

    :param data: The element that you wish to store in the tree
    :param index: The way for you to find this element of data again through
        :func:`Tree.get_node`
    :param level: A way of grouping nodes together, increasing the level and
        then removing all nodes of a given level is a quick way to batch
        remove nodes

    :type data: anything
    :type index: int
    :type level: int

    :attributes:

    `data` : anything
        The element of data being stored in tree
    `index`: int
        The way for you to find this element of data again through
        :func:`Tree.get_node`
    `level`: int
        A way of grouping nodes together, increasing the level and
        then removing all nodes of a given level is a quick way to batch
        remove nodes
    `children`: list of TreeNodes
        Nodes whose parent is this node
    `parent`: TreeNode
        the parent of the current node
    """

    __slots__ = ["data", "index", "level", "children", "parent"]

    def __init__(self, data, index, level):
        super(self.__class__, self).__init__(data, index, level)

    def add_child(self, node):
        """
        adds a new node as a child to the current node

        :param node: the node you want to as a child to this node
        :type node: TreeNodeDynamic
        :return: None
        """

        self.children.append(node)

    def remove_child(self, node):
        """
        removes a node as a child from the current ndoe

        :param node: the node want to remove as a child to this node. must
            already be a child
        :type node: TreeNodeDynamic

        :return:
        """

        if node not in self.children:
            raise exceptions.TreeException(
                "cannot remove child from node as it is not a child")

        self.children.remove(node)


class TreeNodeStatic(TreeNode):
    """
    A TreeNode container that is specific for TreeStatic trees. That is
    it can be connected to an unlimited number of other nodes. There is no
    reason to instantiate this class outside TreeStatic.

    :param data: The element that you wish to store in the tree
    :param index: The way for you to find this element of data again through
        :func:`Tree.get_node`
    :param level: A way of grouping nodes together, increasing the level and
        then removing all nodes of a given level is a quick way to batch
        remove nodes
    :param n_children: The number of connections this node can have, in
        TreeDynamic this does nothing.

    :type data: anything
    :type index: int
    :type level: int
    :type n_children: int

    :attributes:

    `data` : anything
        The element of data being stored in tree
    `index`: int
        The way for you to find this element of data again through
        :func:`Tree.get_node`
    `level`: int
        A way of grouping nodes together, increasing the level and
        then removing all nodes of a given level is a quick way to batch
        remove nodes
    `children`: list of TreeNodes
        Nodes whose parent is this node
    `parent`: TreeNode
        the parent of the current node

    """

    __slots__ = ["data", "index", "level", "children", "parent"]

    def __init__(self, data, index, level, n_children):
        super(self.__class__, self).__init__(data, index, level, n_children)

    def add_child(self, node, pos):
        """
        adds a new node a a child at a specific position in this nodes child's
        list.

        :param node: the node you want to add as a child of this node
        :param pos: the position in this nodes children's list you want to
            place this new child

        :type node: TreeNodeStatic
        :type pos: int

        :return: None
        """

        if pos < 0:
             raise exceptions.TreeEndIndexException(
                self, pos,
                """attempted to add a new child to node: %d at endpos
                    %d but thats lower then zero not allowed"""
                    % (self.index, pos))

        if pos >= len(self.children):
            raise exceptions.TreeEndIndexException(
                self, pos,
                 """attempted to add a new child to node: %d at endpos
                    %d that is larger then the number of children
                    this node has """ % (self.index, pos))

        if self.children[pos] is not None:
            raise exceptions.TreeEndIndexException(
                self, pos,
                """attempted to add a new child to node: %d at endpos
                    %d but it is full""" % (self.index, pos))

        self.children[pos] = node

    def remove_child(self, node):
        """
        Remove a child at a specific position in child list

        :param node: the node to remove
        :type node: TreeNodeStatic

        :return: None
        """

        if node not in self.children:
            raise exceptions.TreeException(
                "cannot remove child from node as it is not a child")

        i = self.children.index(node)
        self.children[i] = None

    def copy(self):
        """
        deep copies tree node and also tries to deep copy data. Will work if
        data element has a function copy().

        :return: a deep copy of node
        :rtype: TreeNodeStatic
        """

        new_data = None
        try:
            new_data = self.data.copy()
        except:
            new_data = self.data
        c = TreeNodeStatic(new_data, self.index, self.level, len(self.children))
        return c


