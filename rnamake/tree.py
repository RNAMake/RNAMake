import Queue

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

    def __init__(self):
        self.nodes = []
        self.last_node = None
        self.index = 0
        self.level = 0

        #iterator stuff
        self.current_node = None

    def get_node(self, index):
        for n in self.nodes:
            if n.index == index:
                return n

        raise ValueError("cannot find node with index "+ str(index))

    def remove_node(self, node=None, index=None):
        if node is None and index is not None:
            node = self.get_node(index)
        elif node is None and index is None:
            raise ValueError("need to include node or index in remove_node")

        if node.parent is not None:
            node.parent.remove_child(node)
        self.nodes.remove(node)
        node.parent = None
        self.last_node = node.parent
        #self.index -= 1

    def remove_node_level(self, level=None):
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
        self.level += 1

    def decrease_level(self):
        self.level -= 1

class TreeDynamic(Tree):
    def __init__(self):
        super(TreeDynamic, self).__init__()

    def add_data(self, data, parent_index=-1):
        n = TreeNodeDynamic(data, self.index, self.level)
        parent = self.last_node
        if parent_index != -1:
            parent = self.get_node(parent_index)

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


