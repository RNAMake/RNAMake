import tree
import Queue
import math
import secondary_structure
from collections import namedtuple

def compare_ss_tree(ss_tree1, ss_tree2):

    bp_score = {
        "G+C-A+U" : 0.25,
        "A+U-G+C" : 0.25,
        "G+C-U+A" : 0.50,
        "U+A-G+C" : 0.50,
        "C+G-A+U" : 0.50,
        "A+U-C+G" : 0.50,
        "C+G-U+A" : 0.25,
        "U+A-C+G" : 0.25
    }

    n_score = {
        "G-A" : 1,
        "G-U" : 2,
        "G-C" : 3,
        "A-G" : 1,
        "A-C" : 2,
        "A-U" : 3,
        "C-U" : 1,
        "C-A" : 2,
        "C-G" : 3,
        "U-C" : 1,
        "U-A" : 3,
        "U-G" : 2
    }

    score = 0

    for n1 in ss_tree1:
        n2 = ss_tree2.get_node(n1.index)

        if n1.data.type == SS_Type.SS_BP:
            if n1.data.seq() ==  n2.data.seq():
                continue
            else:
                str_id = n1.data.seq() + "-" + n2.data.seq()
                if str_id not in bp_score:
                    score += 0.75
                else:
                    score += bp_score[str_id]

        else:
            spl1 = n1.data.seq().split("+")
            spl2 = n2.data.seq().split("+")

            if len(spl1) != len(spl2):
                return 1000

            for i in range(len(spl1)):
                if spl1[i] == spl2[i]:
                    continue

                min = len(spl1[i])
                if len(spl2[i]) < min:
                    min = len(spl2[i])
                score += (math.fabs(len(spl1[i]) - len(spl2[i])))*2

                for j in range(min):
                    if spl1[i][j] == spl2[i][j]:
                        continue
                    c_n_score = n_score[ spl1[i][j] + "-" + spl2[i][j]]
                    score += c_n_score

    return score


class SS_Tree(object):
    def __init__(self, sequence=None, dot_bracket=None, ss=None):
        self.tree = tree.TreeDynamic()

        if sequence is None and dot_bracket is None and ss is None:
            raise ValueError("must supply either sequence and dot_bracket or " +\
                             "secondary structure object")

        if   sequence is not None and dot_bracket is not None:
            self.sequence, self.dot_bracket = sequence, dot_bracket
        elif ss is not None:
            self.sequence, self.dot_bracket = ss.sequence(), ss.dot_bracket()
        else:
            raise ValueError("must supply either sequence and dot_bracket or " +\
                             "secondary structure object")

        if len(self.dot_bracket) != len(self.sequence):
            raise ValueError("sequence and dot bracket are not the same length")

        if self.dot_bracket[0] != '(' and self.dot_bracket[0] != '.':
            raise ValueError("secondary structure is not valid did you flip seq and ss?")

        self._build_tree()

    def get_node(self, i):
        return self.tree.get_node(i)

    def __len__(self):
        return len(self.tree)

    def __iter__(self):
        self.tree.__iter__()
        return self

    def next(self):
        return self.tree.next()

    def seq_from_nodes(self, nodes):
        seen = {}
        for n in nodes:
            for ss_d in n.data.ss_data:
                for i in range(ss_d.bounds[0], ss_d.bounds[1]+1):
                    seen[i] = 1

        start, size, found = -1, 0, 0
        ss_data = []

        for i in range(len(self.seq)):
            if i in seen:
                if found == 0:
                    found, start, size = 1, i, 1
                else:
                    size += 1
            elif(start != -1 and found == 1):
                found = 0
                seq = self.seq[start:start+size]
                ss  = self.ss [start:start+size]
                ss_data.append(SS_Data(seq, ss, [start,start+size-1]))

        if found == 1:
            seq = self.seq[start:start+size]
            ss  = self.ss [start:start+size]
            ss_data.append(SS_Data(seq, ss, [start,start+size-1]))

        return ss_data

    def _build_tree(self):
        xb, yb = 0, len(self.sequence)-1
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
                    ss = current_pair.node.ss
                    bounds = current_pair.node.bounds

                    for n in part_of_nway:
                        for i in range(len(n.ss.chains)):
                            if len(n.ss.chains[i].sequence) > 0:
                                ss.chains.append(n.ss.chains[i])
                                bounds.append(n.bounds[i])

                        if n.type == SS_Type.SS_BULGE:
                            nxb = n.bound_side(0, Bound_Side.RIGHT)+1
                            nyb = n.bound_side(1, Bound_Side.LEFT) -1
                            next_level_2 = self._build_tree_level(nxb, nyb)
                            for n2 in next_level_2:
                                next_level.append(n2)

                    new_node = SS_NodeData(SS_Type.SS_NWAY, ss, bounds)
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
        bounds = [[xb-1,xb-1], [yb+1,yb+1]]
        ss_chains = [secondary_structure.SecondaryStructureChain(),
                     secondary_structure.SecondaryStructureChain()]
        hairpin = 0
        current = None

        if self.dot_bracket[xb] == '.':
            end_x = self._get_dot_bounds(xb, 0)
            for i in range(xb, end_x+1):
                ss_chains[0].sequence    += self.sequence[i]
                ss_chains[0].dot_bracket += self.dot_bracket[i]
            if end_x == yb:
                hairpin = 1
            bounds[0] = [xb, end_x]

        if self.dot_bracket[yb] == '.' and hairpin == 0:
            end_y = self._get_dot_bounds(yb, 1)
            for i in range(end_y, yb+1):
                ss_chains[1].sequence    += self.sequence[i]
                ss_chains[1].dot_bracket += self.dot_bracket[i]
            bounds[1] = [end_y, yb]

        if self.dot_bracket[xb] =='&' or self.dot_bracket[xb] == '+':
            bounds[0] = [xb, xb]
            ss_chains[0].sequences   = '&'
            ss_chains[0].dot_bracket = '&'
            return SS_NodeData(SS_Type.SS_SEQ_BREAK,
                               secondary_structure.SecondaryStructure(ss_chains),
                               bounds)

        type = None
        if len(ss_chains[0].sequence) == 0 and len(ss_chains[1].sequence) == 0:
            pair = self._get_bracket_pair(xb)
            ss_chains[0].sequence = self.sequence[xb]
            ss_chains[0].dot_bracket = "("
            ss_chains[1].sequence = self.sequence[pair]
            ss_chains[1].dot_bracket = ")"
            bounds = [[xb, xb], [pair, pair]]
            type = SS_Type.SS_BP
            if(self.dot_bracket[xb] == '['):
               type = SS_Type.SS_PSEUDO_BP

        elif hairpin:
            type = SS_Type.SS_HAIRPIN

        else:
            type = SS_Type.SS_BULGE

        current = SS_NodeData(type,
                              secondary_structure.SecondaryStructure(ss_chains),
                              bounds)

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
        for i, e in enumerate(self.dot_bracket):
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
            for i in range(pos+1,len(self.dot_bracket)):
                if self.dot_bracket[i] != ".":
                    return i-1
            return len(self.dot_bracket)-1
        else:
            for i in range(pos-1,0,-1):
                if self.dot_bracket[i] != ".":
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
    def __init__(self, type, ss, bounds):
        self.ss, self.type, self.bounds = ss, type, bounds

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

    def sequence(self):
        return self.ss.sequence()

