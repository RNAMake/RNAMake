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
import numpy as np


class NameElements(object):
    def __init__(self, motif_name, helix_direction, start_helix_count,
                 start_index, end_helix_count, end_index, flip_direction):
        self.motif_name, self.helix_direction, self.start_helix_count = \
            motif_name, int(helix_direction), int(start_helix_count)
        self.start_index, self.end_helix_count, self.end_index, self.flip_direction = \
            int(start_index), int(end_helix_count), int(end_index), int(flip_direction)


class MotifTreeState(object):
    def __init__(self, name, start_index, size, score, beads, end, end_index,
                 flip, build_string):
        self.name, self.start_index, self.size = name, start_index, size
        self.score, self.beads, self.end_index = score, beads, end_index
        self.flip, self.build_string, self.end_state = flip, build_string, end


class MotifTreeStateContainer(object):
    def __init__(self,mts,lib_type):
        self.mts = mts
        self.lib_type = lib_type


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
            motif_tree_state = MotifTreeState(name, name_elements.start_index,
                                              size, score, beads, end_state,
                                              name_elements.end_index, flip,
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
                   motif_type.type_to_str(mtype)
            return mtype, path
        else:
            return motif_type.UNKNOWN, libpath


class MotifTreeStateNode(object):
    def __init__(self, mts, level, parent, lib_type, children_lib_types):
        self.mts, self.level, self.parent = mts, level, parent
        self.lib_type, self.children_lib_types = lib_type, children_lib_types
        self.beads, self.score, self.size, self.ss_score = mts.beads, 1000, 0, 0
        self.state = mts.end_state.copy()

    def copy(self):
        c = MotifTreeStateNode(self.mts, self.level, self.parent, self.lib_type,
                               self.children_lib_types)
        c.beads = np.copy(self.beads)
        c.score, c.size, c.ss_score = self.score, self.size, c.ss_score
        c.state = self.state.copy()
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


class MotifTreeStateNodeAligner(object):
    def __init__(self):
        self.r, self.t, self.ref_bp_state = None, None, basepair.ref_bp_state()

    def transform_state(self,parent,child):
        self.r,self.t = parent.state.get_transforming_r_and_t_w_state(self.ref_bp_state)
        self.t += parent.state.d

        child.state = child.mts.end_state.copy()
        new_r,new_d,new_sug = child.state.get_transformed_state(self.r,self.t)
        child.state.set(new_r,new_d,new_sug)
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
        self.aligner = MotifTreeStateNodeAligner()
        self.clash_radius = settings.CLASH_RADIUS
        self.nodes = [ head ]
        self.last_node = head

    def setup_options_and_constraints(self):
        options = { 'sterics'       : 1
                  }
        self.options = option.Options(options)
        self.constraints = {}

    def add_state(self, mts, parent=None):
        if parent is None:
            parent = self.last_node
        new_node = MotifTreeStateNode(mts, parent.level+1, parent, 1, [1])
        self.aligner.transform_state(parent, new_node)
        self.aligner.transform_beads(new_node)
        if self.option('sterics'):
            if new_node.steric_clash():
                return None
        self.nodes.append(new_node)
        self.last_node = new_node
        return new_node

    def to_motiftree(self):
        for i, n in enumerate(self.nodes):
            if i == 0:
                if n.mts.name == "start":
                    mt = motif_tree.MotifTree()
                else:
                    m = motif.str_to_motif(n.mts.build_string)
                    mt = motif_tree.MotifTree(m)
                continue

            m =  motif.str_to_motif(n.mts.build_string)
            parent = mt.nodes [ self.nodes.index(n.parent) ]
            parent_index = n.parent.mts.end_index
            mt_node = mt.add_motif(m, parent=parent, end_index=n.mts.start_index,
                                   end_flip=n.mts.flip, parent_index=parent_index)
            if mt_node is None:
                raise ValueError("could not successfully convert to motiftree")

        return mt

    def _get_default_head(self):
        ref_state = basepair.ref_bp_state()
        start_mts = MotifTreeState("start", 0, 0, 0, [], ref_state, 0, 0, "")
        start_node = MotifTreeStateNode(start_mts, 0, None, 0, [0])
        start_node.beads = []
        return start_node

    def _get_head_node(self, mts):
        start_node = MotifTreeStateNode(mts, 0, None, 0, [0])
        return start_node




def parse_db_name(name):
    spl = name.split("-")
    name_elements = NameElements(*spl)
    return NameElements(*spl)
