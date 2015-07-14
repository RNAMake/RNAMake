import tree
import Queue
import math
import secondary_structure
import uuid
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
            self.ss = secondary_structure.SecondaryStructure(sequence, dot_bracket)
        elif ss is not None:
            self.ss = ss
        else:
            raise ValueError("must supply either sequence and dot_bracket or " +\
                             "secondary structure object")

        self.residues = self.ss.residues()
        pad_residue1 = secondary_structure.Residue("N", ".", 0, "N", uuid.uuid1())
        self.residues.insert(0, pad_residue1)
        pad_residue2 = secondary_structure.Residue("N", ".", len(self.residues),
                                                   "N", uuid.uuid1())
        self.residues.append(pad_residue2)
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

    def _check_from_chain_ends(self, xb, yb):
        if self._is_res_end_of_chain(self.residues[xb].num) or \
           self._is_res_end_of_chain(self.residues[yb].num):
            ss_chains = [secondary_structure.Chain(),
                         secondary_structure.Chain()]
            fake_res1 = secondary_structure.Residue('&', '&', self.residues[xb].num,
                                                    "N", uuid.uuid1())
            fake_res2 = secondary_structure.Residue('&', '&', self.residues[yb].num,
                                                    "N", uuid.uuid1())
            ss_chains[0].residues.append(fake_res1)
            ss_chains[1].residues.append(fake_res2)
            current = SS_NodeData(SS_Type.SS_SEQ_BREAK, ss_chains,
                                  [ self.residues[xb-1], self.residues[yb+1]])
            return current
        return None

    def _build_tree(self):
        xb, yb = 1, len(self.residues)-2
        ss_chains = [secondary_structure.Chain(),
                    secondary_structure.Chain()]
        current = SS_NodeData(SS_Type.SS_SEQ_BREAK, ss_chains,
                              [ self.residues[xb-1], self.residues[yb+1]])
        NodeandIndex = namedtuple('NodeandIndex', 'node index')
        open_nodes = Queue.Queue()
        open_nodes.put(NodeandIndex(current, -1))
        while not open_nodes.empty():
            current_pair = open_nodes.get()

            if current_pair.node.type == SS_Type.SS_HAIRPIN:
                index = self.tree.add_data(current_pair.node, current_pair.index)
                continue

            xb =  self._map_back_to_index(current_pair.node.bound_side(0, Bound_Side.RIGHT))+1
            yb =  self._map_back_to_index(current_pair.node.bound_side(1, Bound_Side.LEFT)) -1


            if xb > yb:
                index = self.tree.add_data(current_pair.node, current_pair.index)
                seq_break = self._check_from_chain_ends(xb-1, yb+1)
                if seq_break is not None:
                    self.tree.add_data(seq_break, index)

                continue

            next_level = self._build_tree_level(xb, yb)

            if len(next_level) == 0:
                seq_break = self._check_from_chain_ends(xb-1, yb+1)
                if seq_break is not None:
                    next_level.append(seq_break)

            if current_pair.node.type == SS_Type.SS_BULGE:
                part_of_nway, not_part_of_nway = [], []
                for n in next_level:
                    if n.type == SS_Type.SS_BULGE or n.type == SS_Type.SS_HAIRPIN:
                        part_of_nway.append(n)
                    else:
                        not_part_of_nway.append(n)
                next_level = not_part_of_nway

                if len(part_of_nway) > 0:
                    ss_chains = current_pair.node.ss_chains
                    bounds = current_pair.node.bounds

                    for n in part_of_nway:
                        for i in range(len(n.ss_chains)):
                            if len(n.ss_chains[i].sequence()) > 0:
                                ss_chains.append(n.ss_chains[i])

                        if n.type == SS_Type.SS_BULGE:
                            nxb = self._map_back_to_index(n.bound_side(0, Bound_Side.RIGHT))+1
                            nyb = self._map_back_to_index(n.bound_side(1, Bound_Side.LEFT)) -1
                            next_level_2 = self._build_tree_level(nxb, nyb)
                            for n2 in next_level_2:
                                next_level.append(n2)

                    new_node = SS_NodeData(SS_Type.SS_NWAY, ss_chains, bounds)
                    index = current_pair.index
                    current_pair = NodeandIndex(new_node, index)

            index = self.tree.add_data(current_pair.node, current_pair.index)
            for n in next_level:
                open_nodes.put(NodeandIndex(n, index))

    def _build_tree_sub(self, xb, yb):
        pass

    def _build_tree_level(self, xb, yb):
        next_level = []

        while xb <= yb:
            if xb == yb and self.residues[xb].dot_bracket != ".":
                break
            child = self._assign_new_node(xb, yb)
            xb =  self._map_back_to_index(child.bound_side(1, Bound_Side.RIGHT))+1
            next_level.append(child)
        return next_level

    def _assign_new_node(self, xb, yb):
        ss_chains = [secondary_structure.Chain(),
                     secondary_structure.Chain()]
        bounds = [self.residues[xb-1], self.residues[yb+1]]
        hairpin = 0
        current = None

        if self.residues[xb].dot_bracket == '.':
            end_x = self._get_dot_bounds(xb, 0)
            for i in range(xb, end_x+1):
                ss_chains[0].residues.append(self.residues[i])
            if end_x == yb:
                hairpin = 1

        if self.residues[yb].dot_bracket == '.' and hairpin == 0:
            end_y = self._get_dot_bounds(yb, 1)
            for i in range(end_y, yb+1):
                ss_chains[1].residues.append(self.residues[i])

        type = None
        if len(ss_chains[0].residues) == 0 and len(ss_chains[1].residues) == 0:
            pair = self._get_bracket_pair(xb)
            ss_chains[0].residues.append(self.residues[xb])
            ss_chains[1].residues.append(self.residues[pair])
            type = SS_Type.SS_BP
            if(self.residues[xb].dot_bracket == '['):
               type = SS_Type.SS_PSEUDO_BP

        elif hairpin:
            type = SS_Type.SS_HAIRPIN

        else:
            type = SS_Type.SS_BULGE

        current = SS_NodeData(type, ss_chains, bounds)

        return current

    def _get_bracket_pair(self, num):
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
        for i, r in enumerate(self.residues):
            if i < num:
                continue
            if r.dot_bracket == "(":
                bracket_count += 1
            if r.dot_bracket == ")":
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

        incr = 1
        if reverse == 1:
            incr = -1

        while 1:
            if pos == 0 or pos == len(self.residues):
                break

            if self.residues[pos].dot_bracket == ".":
                pos += incr
            else:
                pos -= incr
                break



        return pos

    def _is_res_end_of_chain(self, num):
        for c in self.ss.chains:
            if c.last().num == num:
                return 1
        return 0

    def _map_back_to_index(self, r):
        try:
            i = self.residues.index(r)
            return i
        except:
            raise ValueError("cannot find residue " + r )

class SS_Type(object):
   SS_BP        = 0
   SS_BULGE     = 1
   SS_HAIRPIN   = 2
   SS_NWAY      = 3
   SS_PSEUDO_BP = 4
   SS_SEQ_BREAK = 5
   SS_MULTI     = 6

class Bound_Side(object):
    LEFT  = 0
    RIGHT = 1

class SS_NodeData(object):
    def __init__(self, type, ss_chains, bounds):
        self.ss_chains, self.type, self.bounds = ss_chains, type, bounds

    def bound_side(self, pos, side):
        ss_chain = self.ss_chains[pos]
        if len(ss_chain.residues) == 0:
            return self.bounds[pos]

        if side == Bound_Side.LEFT:
            return ss_chain.residues[0]
        else:
            return ss_chain.residues[-1]

    def what(self):
        if   self.type == SS_Type.SS_BP:        return "SS_BP"
        elif self.type == SS_Type.SS_BULGE:     return "SS_BULGE"
        elif self.type == SS_Type.SS_HAIRPIN:   return "SS_HAIRPIN"
        elif self.type == SS_Type.SS_NWAY:      return "SS_NWAY"
        elif self.type == SS_Type.SS_PSEUDO_BP: return "SS_PSEUDO_BP"
        elif self.type == SS_Type.SS_SEQ_BREAK: return "SS_SEQ_BREAK"
        else: raise ValueError("unknown SS_TYPE")

    def sequence(self):
        seqs = [c.sequence() for c in self.ss_chains]
        return "&".join(seqs)

