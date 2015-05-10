import motif_type
import motif_scorer
import settings
import basic_io
import basepair
import base
import option
import settings
import motif_tree
import motif
import util
import residue
import numpy as np
import re


class NameElements(object):
    """
    Stores the results of parsing a databae name for storing a MotifTree

    Attributes
    ----------
    `motif_name` : Str
        The non helix motif stored in the MotifTree
    `helix_direction` : Int

    """
    def __init__(self, motif_name, helix_direction, start_helix_count,
                 start_index, end_helix_count, end_index, flip_direction):
        self.motif_name, self.helix_direction, self.start_helix_count = \
            motif_name, int(helix_direction), int(start_helix_count)
        self.start_index, self.end_helix_count, self.end_index, self.flip_direction = \
            int(start_index), int(end_helix_count), int(end_index), int(flip_direction)

    def get_name(self):
        name = self.motif_name + "-" + str(self.helix_direction) + "-" + \
               str(self.start_helix_count) + "-" + str(self.start_index) + "-" +\
               str(self.end_helix_count) + "-" + str(self.end_index) + "-" +\
               str(self.flip_direction)
        return name


class MotifTreeState(object):
    def __init__(self, name, start_index, size, score, beads, ends, flip, build_string):
        self.name, self.start_index, self.size = name, start_index, size
        self.score, self.beads, self.end_states = score, beads, ends
        self.flip, self.build_string = flip, build_string

    def to_str(self):
        s = self.name + "|" + str(self.start_index) + "|" + str(self.size) + \
            "|" + str(self.score) + "|" + basic_io.points_to_str(self.beads) + \
            "|"
        for state in self.end_states:
            if state is not None:
                s += state.to_str()
            s += "E"
        s += "|" + str(self.flip) + "|" + self.build_string
        return s

    def to_f_str(self):
        s = self.name + "|" + str(self.score) + "|" + str(self.size) + \
            "|" + str(self.flip) +  "|" + self.build_string + "|" + \
             basic_io.points_to_str(self.beads) + "|"

        for end in self.end_states:
            if end is not None:
                s += end.to_str()
            else:
                s += " "
            s += "|"
        s += "\n"
        return s

class MotifTreeStateNodeAligner(object):
    def __init__(self):
        self.r, self.t, self.ref_bp_state = None, None, basepair.ref_bp_state()

    def transform_state(self, parent_end, parent, child):
        self.r,self.t = parent_end.get_transforming_r_and_t_w_state(self.ref_bp_state)
        self.t += parent_end.d

        for i, s in enumerate(child.mts.end_states):
            if s is None:
                continue
            new_r,new_d,new_sug = s.get_transformed_state(self.r,self.t)
            child.states[i].set(new_r,new_d,new_sug)

        child.size = child.mts.size + parent.size
        child.ss_score = child.mts.score + parent.ss_score

    def transform_beads(self,child):
        if len(child.mts.beads) > 0:
            child.beads = np.dot(child.mts.beads, self.r.T) + self.t


def parse_db_name(name):
    spl = name.split("-")
    return NameElements(*spl)


def str_to_motif_tree_state(s):
    spl = s.split("|")
    mts_elements = spl[::]
    mts_elements[1] = int(spl[1])
    mts_elements[2] = float(spl[2])
    mts_elements[3] = float(spl[3])
    mts_elements[4] = basic_io.str_to_points(spl[4])
    states = []
    for bp_str in spl[5].split("E")[:-1]:
        if len(bp_str) < 5:
            states.append(None)
        else:
            states.append(basepair.str_to_basepairstate(bp_str))
    mts_elements[5] = states
    mts_elements[6] = int(spl[6])
    return MotifTreeState(*mts_elements)


def str_to_motif_tree_state_tree(s):
    spl = s.split("#")
    for i, node_str in enumerate(spl[:-1]):
        node_spl = node_str.split("!")
        mts = str_to_motif_tree_state(node_spl[0])
        if i == 0:
            if mts.name == "start":
                mtst = MotifTreeStateTree()
            else:
                mtst = MotifTreeStateTree(mts)
            continue
        parent = mtst.nodes [ int(node_spl[2]) ]
        children_lib_types = [int(x) for x in node_spl[5].split(" ")]
        new_node = MotifTreeStateNode(mts, len(mtst.nodes), parent,
                                      int(node_spl[3]), children_lib_types)
        states = []
        for bp_str in node_spl[4].split("E")[:-1]:
            if len(bp_str) < 5:
                states.append(None)
            else:
                states.append(basepair.str_to_basepairstate(bp_str))
        new_node.beads = basic_io.str_to_points(node_spl[6])
        parent.children[ int(node_spl[1]) ] = new_node
        mtst.nodes.append(new_node)
    return mtst


def ref_mts():
    ref_motif = motif.ref_motif()
    ref_bp = ref_motif.ends[0]
    start_beads = ref_bp.res1.get_beads() + ref_bp.res2.get_beads()
    beads = []
    for b in start_beads:
        if b.btype != residue.BeadType.PHOS:
            beads.append(b.center)
    start_mts = MotifTreeState("start", 1, 0, 0, beads, [ref_bp.state()], 0,
                               ref_motif.to_str())
    return start_mts


def motif_to_state(m, end_index=0, end_flip=0):
    mt = motif_tree.MotifTree()
    mt.add_motif(m, end_index=end_index, end_flip=end_flip)
    m_copy = mt.nodes[1].motif
    name = m_copy.name+"-"+str(end_index)+"-"+str(end_flip)+"-0-0-0-0"
    available_ends = mt.nodes[1].available_ends()
    ends = [ None for e in m_copy.ends]
    end_indexes = []
    for i, end in enumerate(m_copy.ends):
        if end_index == i:
            continue
        ends[i] = end.state()

    beads = []
    for b in m_copy.beads:
        if b.btype != 0:
            beads.append(b.center)

    ms = motif_scorer.MotifScorer()
    score = ms.score(m)
    mts = MotifTreeState(name, end_index, len(m.residues()), score, beads,
                         ends, end_flip, m_copy.to_str())
    return mts


def motif_to_state_simple(m, end_index=0, end_flip=0):
    mt = motif_tree.MotifTree()
    mt.add_motif(m.copy(), end_index=end_index, end_flip=0)
    m_copy = mt.nodes[1].motif
    ends = [ None for e in m_copy.ends ]
    for i, end in enumerate(m_copy.ends):
        if end_index == i:
            continue
        ends[i] = end.state()

    beads = []
    for b in m_copy.beads:
        if b.btype != 0:
            beads.append(b.center)

    ms = motif_scorer.MotifScorer()
    score = ms.score(m)

    mts = MotifTreeState(m_copy.name, end_index,
                         len(m.residues()), score, beads,
                         ends, end_flip, m_copy.to_str())
    return mts


def generate_clash_files(mtype1, mtype2):
    data_path = settings.PRECOMPUTED_PATH + "motif_tree_states/"
    data_path += motif_type.type_to_str(mtype1) + "_"
    data_path += motif_type.type_to_str(mtype2) + ".clist"
    lib1 = MotifTreeStateLibrary(mtype1)
    lib2 = MotifTreeStateLibrary(mtype2)
    mtst = MotifTreeStateTree()

    f = open(data_path,'w')
    mtst = MotifTreeStateTree()
    for mts1 in lib1.motif_tree_states:
        node = mtst.add_state(mts1)
        if node is None:
            continue
        for mts2 in lib2.motif_tree_states:
            node = mtst.add_state(mts2)
            if node is not None:
                mtst.remove_node()
                continue
            f.write(mts1.name + " " + mts2.name + "\n")
        #print len(mtst.nodes)
        mtst.remove_node()

    f.close()


if __name__ == '__main__':
    generate_clash_files(motif_type.TWOWAY, motif_type.TWOWAY)
    #generate_clash_files(motif_type.NWAY, motif_type.TWOWAY)

