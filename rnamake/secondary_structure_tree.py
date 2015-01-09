import sys
import re

def get_bracket_pair(pos,ss):
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
    for i,e in enumerate(ss):
        if e == "(":
            bracket_count += 1
        elif e == ")":
            bracket_count -= 1
            if bracket_count == 0:
                return i
    raise ValueError("cannot find pair")

def get_dot_bounds(pos,ss,reverse=0):
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
        for i in range(pos+1,len(ss)):
            if ss[i] != ".":
                return i-1
        return len(ss)-1
    else:
        for i in range(pos-1,0,-1):
            if ss[i] != ".":
                return i+1
        return 0

    return pos


class SecondaryStructureNode(object):
    """
    Holds an element of secondary structure contained in
    SecondaryStructureTree. This is a base class and is never directly
    instantiated.
    """

    def _assign_children(self):
        """
        assigns the next child node bases on whether there is a "." (Bulge) or
        "(" (Basepair) and passes allong both the secondary structure and
        sequence information for all of its children
        """

        if self.ss[-1] == ".":
            bulge = SSN_Bulge(self.ss,self.seq,parent=self)
            self.children.append(bulge)

            self.ss = []
            self.seq = []

        elif self.ss[0] == "(":

            pair = get_bracket_pair(0,self.ss)
            ##print self.ss[0:pair+1]
            bp = SSN_Basepair(self.ss[0:pair+1],self.seq[0:pair+1],parent=self)
            self.children.append(bp)

            self.ss = self.ss[pair+1:]
            self.seq = self.seq[pair+1:]

        else:

            bulge = SSN_Bulge(self.ss,self.seq,parent=self)
            self.children.append(bulge)

            self.ss = []
            self.seq = []
        return 1

    def assign_all_children(self):
        """
        Each time a secondary structure element is detected,
        either a basepair or a bulge a node is created and the secondary
        structure elements between it are sent to a new node

        Example:

        Current State\n
        self.ss => [(, (, (, (, ., ., ., ), ), ), )]\n
        self.seq => [G, G, G, C, U, U, C, G, G, C, C, C]

        first element in self.ss is "(" so new basepair node is created

        child = Basepair([(, (, (, (, ., ., ., ), ), ), )],[G, G, G, C, U, U, C, G, G, C, C, C])\n

        child.ss => [(, (, (, ., ., ., ), ), )]\n
        child.seq => [G, G, C, U, U, C, G, G, C, C]
        child.bp_type => GC

        Notice that the child.ss and child.seq are two elements shorter. The
        first and last elements were removed from both child.ss and child.seq
        these become the bp type and the rest of the ss and seq are recursively
        shortened as they are passed to this nodes children.
        """

        while len(self.ss) > 0:
            self._assign_children()

    def get_ss_and_seq(self):
        """
        generates both the secondary structure and sequence recursively of
        all children nodes

        :returns: str,str -- secondary structure,sequence
        """
        pass

class SSN_SingleStranded(SecondaryStructureNode):
    """
    Represents single stranded RNA before the first basepair and after the
    last basepair.

    :param seq: Sequence associated with secondary structure
    :type seq: List
    :param parent: parent node
    :type parent: SecondaryStructureNode
    :state parent: None

    """
    def __init__(self,ss,seq,parent=None):
        self.seq = seq
        self.parent = parent
        self.children = []
        self.ss_type = "SingleStranded"

        self.sx = []
        self.sy = []

        self.x_seq = []
        self.y_seq = []

        self.seq = seq
        self.ss = ss

        #print ss,seq

        #remove dots from start
        end = get_dot_bounds(0,self.ss)
        if end > -1 and self.ss[0] == ".":
            for i in range(end+1):
                self.sx.append(self.ss.pop(0))
                self.x_seq.append(self.seq.pop(0))

        end = get_dot_bounds(len(self.ss)-1,self.ss,reverse=1)
        diff = len(self.ss)-end
        #print "".join(self.ss),diff
        if diff > 0 and self.ss[-1] == ".":
            for i in range(diff):
                self.sy.append(self.ss.pop())
                self.y_seq.append(self.seq.pop())

        self.org_x_seq = self.x_seq
        self.org_y_seq = self.y_seq

    def get_ss_and_seq(self):
        ss = ""
        seq = ""
        for c in self.children:
            c_ss,c_seq = c.get_ss_and_seq()
            ss += c_ss
            seq += c_seq

        return  "".join(self.sx) + ss + "".join(self.sy), "".join(self.x_seq) + seq + "".join(self.y_seq[::-1])

    def revert_sequence(self):
        for i,e in enumerate(self.org_x_seq):
            self.x_seq[i] = e
        for i,e in enumerate(self.org_y_seq):
            self.y_seq[i] = e

    def set_seq(self,seq):
        seqs = re.split("\-",seq)
        self.x_seq = list(seqs[0])
        if len(seqs) > 1:
            self.y_seq = list(seqs[1])

class SSN_Basepair(SecondaryStructureNode):
    """
    Holds a basepair

    :param ss: Secondary structure in dot bracket notation
    :type ss: List
    :param seq: Sequence associated with secondary structure
    :type seq: List
    :param parent: parent node
    :type parent: SecondaryStructureNode
    :state parent: None

    Attributes
    ----------
    `bp_type` : str
        the concatenation of both residue strings, ex. AU or GC or CG etc
    `children` :  List of SecondaryStructureNodes
        nodes generated from the ss and seq variables
    `res1` : str
        residue from left side of secondary structure
    `res2` : str
        residue from right side of secondary structure
    `ss` : List
        secondary structure
    `seq` : List
        sequence
    `ss_type` : str
        what type of SSN is this
    """

    def __init__(self,ss,seq,parent=None):
        self.parent = parent
        self.children = []
        self.ss_type = "Basepair"
        self.seq = seq
        self.ss = ss

        self.res1 = self.seq.pop(0)
        self.res2 = self.seq.pop()
        self.bp_type = self.res1+self.res2

        self.ss.pop(0)
        self.ss.pop()

    def get_ss_and_seq(self):
        ss = ""
        seq = ""
        for c in self.children:
            c_ss,c_seq = c.get_ss_and_seq()
            ss += c_ss
            seq += c_seq

        return "(" + ss + ")", self.res1 + seq + self.res2

    def set_seq(self,seq):
        self.res1 = seq[0]
        self.res2 = seq[1]
        self.bp_type = seq

class SSN_Bulge(SecondaryStructureNode):
    """
    Holds either Bulge,Junction or Hairpin probably could be renamed

    :param ss: Secondary structure in dot bracket notation
    :type ss: List
    :param seq: Sequence associated with secondary structure
    :type seq: List
    :param parent: parent node
    :type parent: SecondaryStructureNode
    :state parent: None

    Attributes
    ----------
    `bp_type` : str
        the concatenation of both residue strings, ex. AU or GC or CG etc
    `children` :  List of SecondaryStructureNodes
        nodes generated from the ss and seq variables
    `sx` : List
        secondary structure from left side until next basepair
    `sy` : List
        secondary structure from right side until next basepair
    `x_seq` : List
        sequence from left side until next basepair
    `y_seq` : List
        sequence from right side until next basepair
    `ss` : List
        secondary structure
    `seq` : List
        sequence
    `ss_type` : str
        what type of SSN is this
    """
    def __init__(self,ss,seq,parent=None):
        self.parent = parent
        self.children = []
        self.ss_type = "Bulge"

        self.sx = []
        self.sy = []

        self.x_seq = []
        self.y_seq = []

        self.seq = seq
        self.ss = ss

        #print "start"
        #print "".join(ss)
        #print "".join(seq)

        #remove dots from start
        end = get_dot_bounds(0,self.ss)
        if end > -1 and self.ss[0] == ".":
            for i in range(end+1):
                self.sx.append(self.ss.pop(0))
                self.x_seq.append(self.seq.pop(0))

        end = get_dot_bounds(len(self.ss)-1,self.ss,reverse=1)
        diff = len(self.ss)-end
        #print "".join(self.ss),diff
        if diff > 0 and self.ss[-1] == ".":
            for i in range(diff):
                self.sy.append(self.ss.pop())
                self.y_seq.append(self.seq.pop())

        #print self.x_seq,self.y_seq

        #check for n junction
        #for e in ss:

    def get_ss_and_seq(self):
        ss = ""
        seq = ""
        for c in self.children:
            c_ss,c_seq = c.get_ss_and_seq()
            ss += c_ss
            seq += c_seq

        return  ("".join(self.sx) + ss + "".join(self.sy),
                "".join(self.x_seq) + seq + "".join(self.y_seq[::-1]))

    def set_seq(self,seq):
        seqs = re.split("\-",seq)
        self.x_seq = list(seqs[0])
        if len(seqs) > 1:
            self.y_seq = list(seqs[1])


class SecondaryStructureTree(object):
    """
    Converts a dot bracket notation for secondary structure to a tree
    format. The sequence corresponding to the secondary structure is also
    required to allow for easy evaluations of what are the contents of each
    node.

    :param ss: Secondary structure in dot bracket notation
    :type ss: str
    :param seq: Sequence associated with secondary structure
    :type seq: str
    :raises: ValueError

    Example:
    ::
        ss_tree = SecondaryStructureTree("(.((....))..)", "GUGCUUCGGCUUC")

    Here would be the iternal node structure

    Internal node structure would be:
    ::
        SSN_Basepair G-C
            SSN_Bulge U-UU
                SSN_Basepair G-C
                    SSN_Basepair C-G
                        SSN_Bulge UUCG-


    Attributes
    ----------
    `basepairs` : List of Basepair Objects
        all basepair nodes
    `bulges` : List of Bulge Objects
        all bulge nodes
    `nodes` : List of SS Node Objects
        all nodes
    `ss` : List
        secondary structure
    `seq` : List
        sequence

    """

    def __init__(self,ss=None,seq=None):
        """
        Generates a new SecondaryStructureTree

        :param ss: Secondary structure in dot bracket notation
        :type ss: str
        :param seq: Sequence associated with secondary structure
        :type seq: str
        :raises: ValueError
        """

        self.basepairs = []
        self.bulges = []
        self.single_strands = []
        self.ss = list(ss)
        self.seq = list(seq)
        self.nodes_to_optimize = []
        self.nodes = []

        #empty tree, currently used for making decoy secondary structure /
        #sequences
        if ss == None:
            return

        if len(ss) != len(seq):
            raise ValueError("secondary structure and sequence are not the same length")

        #remove the single stranded start, these parts are NOT included in the
        #node tree
        #self._remove_single_strand_start()

        #starts recursive process of building nodes
        while len(self.ss) > 0:
            self._build_tree()

        for node in self.nodes:
            if node.ss_type == "Basepair":
                self.basepairs.append(node)
            elif node.ss_type == "Bulge":
                self.bulges.append(node)
            else:
                self.single_strands.append(node)

        for node in self.basepairs:
            if node.res1 == "N" or node.res2 == "N":
                self.nodes_to_optimize.append(node)
        for node in self.bulges:
            string = node.x_seq + node.y_seq
            for e in string:
                if e == "N":
                    self.nodes_to_optimize.append(node)
                    break
        for node in self.single_strands:
            string = node.x_seq + node.y_seq
            for e in string:
                if e == "N":
                    self.nodes_to_optimize.append(node)
                    break

        #for node in self.

    def _remove_single_strand_start(self):
        """
        Removes the start and end single stranded portions of the dot bracket
        notation and are put in SingleStranded nodes

        Example:

        self.ss => [., ., ., ., (, (, (, (, ., ., ., ), ), ), ), ., ., .]
        self.seq => [A, A, A, A, G, G, G, C, U, U, C, G, G, C, C, C, A, A, A]

        >>>self._remove_single_strand_start()

        self.ss => [(, (, (, (, ., ., ., ), ), ), )]
        self.seq => [G, G, G, C, U, U, C, G, G, C, C, C]

        """

        if self.ss[0] == ".":
            end = get_dot_bounds(0,self.ss)
            ss = SSN_SingleStranded(self.seq[0:end+2])
            self.nodes.append(ss)
            for i in range(end+1):
                self.ss.pop(0)
                self.seq.pop(0)
            ##print self.ss

        if self.ss[-1] == ".":
            end = get_dot_bounds(len(self.ss)-1,self.ss,reverse=1)
            ss = SSN_SingleStranded(self.seq[end-1:])
            self.nodes.append(ss)
            diff = len(self.ss)-end
            for i in range(diff):
                self.ss.pop()
                self.seq.pop()

    def _build_tree(self):
        """
        Starts the recursive process of building the tree from the secondary
        structure list. Each time a secondary structure element is detected,
        either a basepair or a bulge a node is created and the secondary
        structure elements between it are sent to a new node

        Example:

        Current State:
        self.ss => [(, (, (, (, ., ., ., ), ), ), )]
        self.seq => [G, G, G, C, U, U, C, G, G, C, C, C]

        first element in self.ss is "(" so new basepair node is created

        node = Basepair([(, (, (, (, ., ., ., ), ), ), )],[G, G, G, C, U, U, C, G, G, C, C, C])

        node.ss => [(, (, (, ., ., ., ), ), )]
        node.seq => [G, G, C, U, U, C, G, G, C, C]
        node.bp_type => GC

        Notice that the node.ss and node.seq are two elements shorter. The
        first and last elements were removed from both node.ss and node.seq
        these become the bp type and the rest of the ss and seq are recursively
        shortened as they are passed to this nodes children.
        """

        node = None
        if self.ss[0] == "." or self.ss[-1] == "." and len(self.nodes) == 0:
            node = SSN_SingleStranded(self.ss,self.seq)
            self.ss = []
            self.seq = []

        elif self.ss[0] == "(":
            pair = get_bracket_pair(0,self.ss)
            node = SSN_Basepair(self.ss[0:pair+1],self.seq[0:pair+1])

            self.ss = self.ss[pair+1:]
            self.seq = self.seq[pair+1:]


        elif self.ss[0] == ".":
            node = SSN_Bulge(self.ss,self.seq)
            self.ss = []
            self.seq = []


        #breath first transverse to vist all node children
        open_nodes = [node]
        while len(open_nodes) > 0:
            current = open_nodes.pop()
            if current in self.nodes:
                continue
            self.nodes.append(current)
            current.assign_all_children()

            open_nodes.extend(current.children)

    def get_ss_and_seq(self):
        return self.nodes[0].get_ss_and_seq()


def main():
    ss  = ".....((..((.(......)))....(((((((....).)))..))).(((....).))..)).....(((((((....)))))))....................."
    seq = "GGAAAGCAAGGACGAAUAAGCCAUAACCAGAGCGAAAGACUCAAUGGAGCCGAAAGAGCAAGCAAUAACUGAUGCUUCGGCAUCAGAAAAGAAACAACAACAACAAC"

    ss_tree = SecondaryStructureTree(ss,seq)

    nss,nseq = ss_tree.get_ss_and_seq()

    #print seq
    #print nseq

    #print ss
    #print nss



if __name__ == '__main__':
    main()



