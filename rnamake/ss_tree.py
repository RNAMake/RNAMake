import tree
import Queue
from collections import namedtuple

class SS_Tree(object):
    def __init__(self, ss, seq):
        self.tree = tree.TreeDynamic()
        self.seq, self.ss = seq, ss
        self._build_tree()

    def __len__(self):
        return len(self.tree)

    def __iter__(self):
        self.tree.__iter__()
        return self

    def next(self):
        return self.tree.next()

    def _build_tree(self):
        xb, yb = 0, len(self.seq)-1
        current = self._assign_new_node(xb, yb)
        index = 0
        NodeandIndex = namedtuple('NodeandIndex', 'node index')
        open_nodes = Queue.Queue()
        open_nodes.put(NodeandIndex(current, -1))
        while not open_nodes.empty():
            current_pair = open_nodes.get()

            if current_pair.node.type == SS_Type.SS_HAIRPIN:
                index = self.tree.add_data(current_pair.node, current_pair.index)
                continue

            xb = current_pair.node.bound_side(0, Bound_Side.RIGHT)+1
            yb = current_pair.node.bound_side(1, Bound_Side.LEFT) -1

            next_level = self._build_tree_level(xb, yb)

            if current_pair.node.type == SS_Type.SS_BULGE:
                part_of_nway, not_part_of_nway = [], []
                for n in next_level:
                    if n.type == SS_Type.SS_BULGE or n.type == SS_Type.SS_HAIRPIN:
                        part_of_nway.append(n)
                    else:
                        not_part_of_nway.append(n)
                next_level = not_part_of_nway

                if len(part_of_nway) > 0:
                    seqs = []
                    bounds = []

                    seqs.append(current_pair.node.seqs[0])
                    bounds.append(current_pair.node.bounds[0])

                    seqs.append(current_pair.node.seqs[1])
                    bounds.append(current_pair.node.bounds[1])

                    for n in part_of_nway:
                        child_seq = n.seqs
                        for i in range(len(child_seq)):
                            if len(child_seq[i]) > 0:
                                seqs.append(child_seq[i])
                                bounds.append(n.bounds[i])

                        if n.type == SS_Type.SS_BULGE:
                            nxb = n.bound_side(0, Bound_Side.RIGHT)+1
                            nyb = n.bound_side(1, Bound_Side.LEFT) -1
                            next_level_2 = self._build_tree_level(nxb, nyb)
                            for n2 in next_level_2:
                                next_level.append(n2)

                    new_node = SS_NodeData(seqs, SS_Type.SS_NWAY, bounds)
                    index = current_pair.index
                    current_pair = NodeandIndex(new_node, index)

            index = self.tree.add_data(current_pair.node, current_pair.index)
            for n in next_level:
                open_nodes.put(NodeandIndex(n, index))


    def _build_tree_level(self, xb, yb):
        next_level = []
        while xb <= yb:
            child = self._assign_new_node(xb, yb)
            xb = child.bound_side(1, Bound_Side.RIGHT)+1
            next_level.append(child)

        return next_level

    def _assign_new_node(self, xb, yb):
        seq = ["", ""]
        bounds = [[xb-1,xb-1], [yb+1,yb+1]]
        hairpin = 0
        current = None

        print xb, yb

        if self.ss[xb] == '.':
            end_x = self._get_dot_bounds(xb, 0)
            for i in range(xb, end_x+1):
                seq[0] += self.seq[i]
            if end_x == yb:
                hairpin = 1
            bounds[0] = [xb, end_x]

        if self.ss[yb] == '.' and hairpin == 0:
            end_y = self._get_dot_bounds(yb, 1)
            for i in range(end_y, yb+1):
                seq[1] += self.seq[i]
            bounds[1] = [end_y, yb]

        if self.ss[xb] =='&' or self.ss[xb] == '+':
            bounds[0] = [xb, xb]
            seq = ['+', '']
            return SS_NodeData(seq, SS_Type.SS_SEQ_BREAK, bounds)

        if len(seq[0]) == 0 and len(seq[1]) == 0:
            pair = self._get_bracket_pair(xb)
            seq[0] = self.seq[xb]
            seq[1] = self.seq[pair]
            bounds[0] = [xb, xb]
            bounds[1] = [pair, pair]

            if(self.ss[xb] == '('):
                current = SS_NodeData(seq, SS_Type.SS_BP, bounds)
            else:
                current = SS_NodeData(seq, SS_Type.SS_PSEUDO_BP, bounds)

        elif hairpin:
            current = SS_NodeData(seq, SS_Type.SS_HAIRPIN, bounds)

        else:
            current = SS_NodeData(seq, SS_Type.SS_BULGE, bounds)

        return current

    def _get_bracket_pair(self, pos):
        """
        finds the position of ")" pair to the current "(" in the secondary
        structure array

        :params pos: position in secondary structure of current "("
        :type pos: int
        :params ss: secondary structure
        :type ss: List
        :returns: int -- the position of the paired ")" in list ss
        :raises: ValueError

        >>>print get_bracket_pair(0, [(,(, ., ), )])
        4
        """
        bracket_count = 0
        for i, e in enumerate(self.ss):
            if i < pos:
                continue
            if e == "(":
                bracket_count += 1
            elif e == ")":
                bracket_count -= 1
                if bracket_count == 0:
                    return i

        raise ValueError("cannot find pair")

    def _get_dot_bounds(self, pos,reverse=0):
        """
        finds the bounds of a current dot stretch

        :params pos: position in secondary structure of current "."
        :type pos: int
        :params ss: secondary structure
        :type ss: List
        :params reverse: direction of search::
            0 -- postive, going left to right in the array
            1 -- negative, going right to left in the array
        :type reverse: int
        :state reverse: 0
        :returns: int -- the position of last "." in stretch

        >>>print get_dot_bounds(1,[(, ., ., ., )])
        3

        >>>print get_dot_bounds(3,[(, ., ., ., )],reverse=1)
        1
        """

        if not reverse:
            for i in range(pos+1,len(self.ss)):
                if self.ss[i] != ".":
                    return i-1
            return len(self.ss)-1
        else:
            for i in range(pos-1,0,-1):
                if self.ss[i] != ".":
                    return i+1
            return 0

        return pos



class SS_Type(object):
   SS_BP        = 0
   SS_BULGE     = 1
   SS_HAIRPIN   = 2
   SS_NWAY      = 3
   SS_PSEUDO_BP = 4
   SS_SEQ_BREAK = 5

class Bound_Side(object):
    LEFT  = 0
    RIGHT = 1

class SS_NodeData(object):
    def __init__(self, seqs, type, bounds):
        self.seqs, self.type, self.bounds = seqs, type, bounds

    def bound_side(self, pos, side):
        bound = self.bounds[pos]
        return bound[side]

    def what(self):
        if   self.type == SS_Type.SS_BP:        return "SS_BP"
        elif self.type == SS_Type.SS_BULGE:     return "SS_BULGE"
        elif self.type == SS_Type.SS_HAIRPIN:   return "SS_HAIRPIN"
        elif self.type == SS_Type.SS_NWAY:      return "SS_NWAY"
        elif self.type == SS_Type.SS_PSEUDO_BP: return "SS_PSEUDO_BP"
        elif self.type == SS_Type.SS_SEQ_BREAK: return "SS_SEQ_BREAK"
        else: raise ValueError("unknown SS_TYPE")


