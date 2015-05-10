import motif_tree_state
import motif_type
import basic_io
import settings
import basepair

class MotifTreeStateLibrary(object):
    def __init__(self, mtype=None, libpath=None, exclude=[], new=0):
        mtype, path = self._parse_args(mtype, libpath)
        self.mtype, self.neighbor_libs, self.children = mtype, [], []
        self.clashes, self.index = {}, 0
        self.new = new
        self.motif_tree_states_dict = self._load_states_from_file(path)
        """motif_tree_states = []
        for mts in self.motif_tree_states:
            s = str(mts.start_index)+str(mts.flip)
            if s in exclude:
                continue
            motif_tree_states.append(mts)
        self.motif_tree_states = motif_tree_states"""

    def motif_tree_states(self):
        return self.motif_tree_states_dict.values()

    def get_state(self, name):
        if name in self.motif_tree_states_dict:
            return self.motif_tree_states_dict[name]
        else:
            return None

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

        motif_tree_state_dict= {}
        for l in lines:
            spl = l.split("|")
            name, score, size = spl[0], float(spl[1]), float(spl[2])
            flip, build_string = int(spl[3]), spl[4]
            beads = basic_io.str_to_points(spl[5])
            end_states = []
            for i in range(6, len(spl)-1):
                if len(spl[i]) < 5:
                    end_states.append(None)
                else:
                    end_states.append(basepair.str_to_basepairstate(spl[i]))
            if not self.new:
                name_elements = motif_tree_state.parse_db_name(name)
                mts = motif_tree_state.MotifTreeState(name,
                                                name_elements.start_index,
                                                size, score, beads, end_states,
                                                flip, build_string)
            else:
                spl = name.split("-")
                start_index = int(spl[1])
                mts = motif_tree_state.MotifTreeState(name,
                                                start_index, size, score, beads,
                                                end_states, flip, build_string)


            motif_tree_state_dict[mts.name] = mts

        return motif_tree_state_dict

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
