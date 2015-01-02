import motif_type
import settings
import basic_io
import basepair


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


def parse_db_name(name):
    spl = name.split("-")
    name_elements = NameElements(*spl)
    return NameElements(*spl)
