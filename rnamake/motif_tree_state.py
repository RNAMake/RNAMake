import motif_type
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


class NameElements(object):
    def __init__(self, motif_name, helix_direction, start_helix_count,
                 start_index, end_helix_count, end_index, flip_direction):
        self.motif_name, self.helix_direction, self.start_helix_count = \
            motif_name, int(helix_direction), int(start_helix_count)
        self.start_index, self.end_helix_count, self.end_index, self.flip_direction = \
            int(start_index), int(end_helix_count), int(end_index), int(flip_direction)


class MotifTreeState(object):
    def __init__(self, name, start_index, size, score, beads, ends, end_indexes,
                 flip, build_string):
        self.name, self.start_index, self.size = name, start_index, size
        self.score, self.beads, self.end_indexes = score, beads, end_indexes
        self.flip, self.build_string, self.end_states = flip, build_string, ends

    def to_str(self):
        s = self.name + "|" + str(self.start_index) + "|" + str(self.size) + \
            "|" + str(self.score) + "|" + basic_io.points_to_str(self.beads) + \
            "|"
        for state in self.end_states:
            s += state.to_str() + "E"
        s += "|"
        for i in self.end_indexes:
            s += str(i) + "E"
        s += "|" + str(self.flip) + "|" + self.build_string
        return s


class MotifTreeStateLibrary(object):
    def __init__(self, mtype=None, libpath=None):
        mtype, path = self._parse_args(mtype, libpath)
        self.mtype, self.neighbor_libs, self.children = mtype, [], []
        self.clashes, self.index = {}, 0
        self.motif_tree_states = self._load_states_from_file(path)

    def get_state(self, name):
        for mts in self.motif_tree_states:
            if mts.name == name:
                return mts
        return None

    def possible_children(self, current_mts):
        self.children = []
        clash_key = ""
        for nlib in self.neighbor_libs:
            for mts in nlib.motif_tree_states:
                clash_key = current_mts.name + " " + mts.name
                if clash_key in self.clashes:
                    continue
                self.children.append(MotifTreeStateContainer(mts, nlib.index))

    def add_neighbor(self, neighbor):
        self.neighbor_libs.append(neighbor)
        path = settings.RESOURCES_PATH + "/precomputed/motif_tree_states/"
        path += motif_type.type_to_str(self.mtype) + "_" + \
                motif_type.type_to_str(neighbor.mtype)
        self.add_clash_file(path)

    def add_clash_file(self, cfile):
        try:
            f = open(cfile)
            lines = f.readlines()
            f.close()
        except IOError:
            raise IOError("cannot open clash file " + cfile)

        for l in lines:
            self.clashes[l.rstrip()] = 1

    def _load_states_from_file(self, file_path):
        f = open(file_path)
        lines = f.readlines()
        f.close()

        motif_tree_states = []
        for l in lines:
            spl = l.split("|")
            name, score, size = spl[0], float(spl[1]), float(spl[2])
            flip, build_string = int(spl[3]), spl[4]
            beads = basic_io.str_to_points(spl[5])
            end_state = basepair.str_to_basepairstate(spl[6])
            name_elements = parse_db_name(name)
            # TODO go back and make sure this still works with multiple ends
            motif_tree_state = MotifTreeState(name, name_elements.start_index,
                                              size, score, beads, [end_state],
                                              [name_elements.end_index], flip,
                                              build_string)
            motif_tree_states.append(motif_tree_state)

            if len(spl[7]) > 1:
                raise ValueError("not implemented yet")

        return motif_tree_states

    def _parse_args(self, mtype, libpath):
        if   mtype is None and libpath is None:
            raise ValueError("must supply motif type or path to the library")
        elif mtype is not None and libpath is not None:
            raise ValueError("cannot supply both mtype and libpath")
        elif mtype is not None:
            motif_type.is_valid_motiftype(mtype)
            path = settings.RESOURCES_PATH + "/precomputed/motif_tree_states/" +\
                   motif_type.type_to_str(mtype) + ".new.me"
            return mtype, path
        else:
            return motif_type.UNKNOWN, libpath


class MotifTreeStateNode(object):
    def __init__(self, mts, index, level, lib_type, children_lib_types):
        self.mts, self.level = mts, level
        self.lib_type, self.children_lib_types = lib_type, children_lib_types
        self.beads, self.score, self.size, self.ss_score = mts.beads, 1000, 0, 0
        self.index = index
        self.states = [ s.copy() for s in self.mts.end_states ]
        self.connections = []
        self.end_status = [1 for index in mts.end_indexes]

    def copy(self):
        c = MotifTreeStateNode(self.mts, self.index, self.level, self.lib_type,
                               self.children_lib_types)
        c.beads = np.copy(self.beads)
        c.score, c.size, c.ss_score = self.score, self.size, c.ss_score
        c.states = [s.copy() for s in self.states]
        return c

    def steric_clash(self):
        dist = 0
        current = self.parent
        while current != None:
            for b1 in self.beads:
                for b2 in current.beads:
                    dist = util.distance(b1, b2)
                    if dist < settings.CLASH_RADIUS:
                        return 1
            current = current.parent
        return 0

    def available_ends(self):
        states = []
        for i, state in enumerate(self.states):
            if self.end_status[i] == 1:
                states.append(state)
        return states

    def parent(self):
        for c in self.connections:
            if c.child == self:
                return c.parent

    def parent_end_index(self):
        for c in self.connections:
            if c.child == self:
                i = c.parent.states.index(c.parent_end)
                return c.parent.mts.end_indexes[i]

    def parent_end(self):
        for c in self.connections:
            if c.child == self:
                return c.parent.states.index(c.parent_end)

    def to_str(self):
        parent_index = -1
        if self.parent() is not None:
            parent_index = self.parent().index
        s = self.mts.to_str() + "!" + str(self.parent_end_index()) + "!" + \
            str(parent_index) + "!" + str(self.lib_type) + "!"
        for state in self.states:
             s += state.to_str() + "E"
        s +="!" + basic_io.point_to_str(self.children_lib_types) + "!"
        s +=basic_io.points_to_str(self.beads)
        return s


class MotifTreeStateConnection(object):
    def __init__(self, parent, child, parent_end):
        self.parent, self.child, self.parent_end = parent, child, parent_end
        pe_index = self.parent.states.index(parent_end)
        self.parent.end_status[pe_index] = 0
        self.parent.connections.append(self)
        self.child.connections.append(self)

    def disconnect(self):
        pe_index = self.parent.states.index(parent_end)
        self.parent.end_status[self.parent_end] = 1
        self.parent.connections.remove(self)
        self.child.connections.remove(self)


class MotifTreeStateNodeAligner(object):
    def __init__(self):
        self.r, self.t, self.ref_bp_state = None, None, basepair.ref_bp_state()

    def transform_state(self, parent_end, parent, child):
        self.r,self.t = parent_end.get_transforming_r_and_t_w_state(self.ref_bp_state)
        self.t += parent_end.d

        for i, s in enumerate(child.mts.end_states):
            new_r,new_d,new_sug = s.get_transformed_state(self.r,self.t)
            child.states[i].set(new_r,new_d,new_sug)

        child.size = child.mts.size + parent.size
        child.ss_score = child.mts.score + parent.ss_score

    def transform_beads(self,child):
        child.beads = np.dot(child.mts.beads, self.r.T) + self.t


class MotifTreeStateTree(base.Base):
    def __init__(self, head_state=None, **options):
        if head_state is None:
            head = self._get_default_head()
        else:
            head = self._get_head_node(head_state)

        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.aligner = MotifTreeStateNodeAligner()
        self.clash_radius = settings.CLASH_RADIUS
        self.nodes = [ head ]
        self.last_node = head

    def setup_options_and_constraints(self):
        options = { 'sterics'       : 1
                  }
        self.options = option.Options(options)
        self.constraints = {}

    def add_state(self, mts, parent=None, parent_end=None):
        if parent is None:
            parent = self.last_node
        if parent_end is not None:
            parent_ends = [ parent_end ]
        else:
            parent_ends = parent.available_ends()

        new_node = MotifTreeStateNode(mts, len(self.nodes), parent.level+1, 1, [1])
        success = 0
        for pe in parent_ends:
            self.aligner.transform_state(pe, parent, new_node)
            self.aligner.transform_beads(new_node)
            if self.option('sterics') == 1:
                if self._steric_clash(new_node):
                    continue
            MotifTreeStateConnection(parent, new_node, pe)
            success=1
            break
        if not success:
            return None
        self.nodes.append(new_node)
        self.last_node = new_node
        return new_node

    def to_motiftree(self, **options):
        for i, n in enumerate(self.nodes):
            if i == 0:
                if n.mts.name == "start":
                    mt = motif_tree.MotifTree(**options)
                else:
                    m = motif.str_to_motif(n.mts.build_string)
                    mt = motif_tree.MotifTree(m, **options)
                continue

            m =  motif.str_to_motif(n.mts.build_string)

            parent = mt.nodes [ self.nodes.index(n.parent()) ]
            parent_index = n.parent_end_index()
            # print parent, n.parent().mts.start_index, parent_index, n.mts.start_index
            mt_node = mt.add_motif(m, parent=parent, end_index=n.mts.start_index,
                                   end_flip=n.mts.flip, parent_index=parent_index)
            if mt_node is None:
                print i, n.mts.name
                raise ValueError("could not successfully convert to motiftree")

        return mt

    def _get_default_head(self):
        start_mts = ref_mts()
        start_node = MotifTreeStateNode(start_mts, 0, 0, 0, [0])
        start_node.index = 0
        return start_node

    def _steric_clash(self, new_node):
        dist = 0
        for n in self.nodes[::-1]:
            for b1 in new_node.beads:
                for b2 in n.beads:
                    dist = util.distance(b1, b2)
                    if dist < self.clash_radius:
                        return 1
        return 0

    def _get_head_node(self, mts):
        start_node = MotifTreeStateNode(mts, 0, 0, 0, [0])
        start_node.index = 0
        return start_node

    def to_str(self):
        s = ""
        for n in self.nodes:
            s += n.to_str() + "#"
        return s

    def to_pdb(self):
        mt = self.to_motiftree()
        mt.to_pdb("mtst.pdb")

    def nodes_to_pdbs(self):
        mt = self.to_motiftree()
        mt.write_pdbs()


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
    states = [ basepair.str_to_basepairstate(bp_str) for bp_str in spl[5].split("E")[:-1]]
    mts_elements[5] = states
    mts_elements[6] = [int(x) for x in spl[6].split("E")[:-1]]
    mts_elements[7] = int(spl[7])
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
        print int(node_spl[2])
        parent = mtst.nodes [ int(node_spl[2]) ]
        children_lib_types = [int(x) for x in node_spl[5].split(" ")]
        new_node = MotifTreeStateNode(mts, len(mtst.nodes), parent.level+1,
                                      int(node_spl[3]), children_lib_types)
        new_node.states = [ basepair.str_to_basepairstate(s) for s in node_spl[4].split("E")[:-1]]
        new_node.beads = basic_io.str_to_points(node_spl[6])
        parent_end = parent.states[ int(node_spl[1]) ]
        MotifTreeStateConnection(parent, new_node, parent_end)
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
    start_mts = MotifTreeState("start", 0, 0, 0, beads, [ref_bp.state()], [0], 0,
                               ref_motif.to_str())
    return start_mts


def motif_to_state(m, end_index=0, end_flip=0):
    mt = motif_tree.MotifTree()
    mt.add_motif(m, end_index=end_index, end_flip=end_flip)
    m_copy = mt.nodes[1].motif
    name = n.name+"-"+str(end_index)+"-"+str(end_flip)
    available_ends = mt.nodes[1].available_ends()
    ends = []
    end_indexes = []
    for e in available_ends:
        ends.append(e.state())
        index = m_copy.ends.index(e)
        end_indexes.append(index)

    mts = MotifTreeState(name, end_index, len(m.residues()), 0, m_copy.beads,
                         ends, end_indexes, end_flip, m_copy.to_str())
    pass
