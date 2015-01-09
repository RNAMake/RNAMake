import rnamake.motif_tree_state as motif_tree_state
import rnamake.settings as settings
import os

class MotifState(object):
    def __init__(self, mts, population):
        self.mts, self.population = mts, population

class MotifEnsemble(object):
    def __init__(self, lib_path=None, start_index=None, flip_direction=None,
                 limit=9999):
        if lib_path == None:
            self.motif_states = []
            return
        self._setup(lib_path, start_index, flip_direction)

    def _get_populations(self, pop_path):
        try:
            f = open(pop_path)
            lines = f.readlines()
            f.close()
        except IOError:
            raise IOError("cannot open MotifStateDistribution population file")

        pops = {}
        for l in lines:
            spl = l.split()
            pops[ spl[0] ] = float(spl[-1])
        return pops

    def _setup(self, lib_path, start_index, flip_direction):
        prediction_dir = settings.RESOURCES_PATH + "/prediction/"
        if os.path.isfile(prediction_dir + "/" + lib_path + ".new.me"):
            mts_lib_path = prediction_dir + "/" + lib_path + ".new.me"
            pop_path = prediction_dir + "/" + lib_path + ".pop"
        elif os.path.isfile(lib_path + ".new.me"):
            raise ValueError("not implemented")

        pops = self._get_populations(pop_path)
        mts_lib = motif_tree_state.MotifTreeStateLibrary(libpath=mts_lib_path)
        self.motif_states = []
        for mts in mts_lib.motif_tree_states:
            if start_index is not None and start_index != mts.start_index:
                continue
            if flip_direction is not None and flip_direction != mts.flip:
                continue
            name_elements = motif_tree_state.parse_db_name(mts.name)
            mname = name_elements.motif_name
            try:
                m_pop = pops[mname]
            except:
                raise ValueError("cannot find "+mname+" in populations")
            motif_state = MotifState(mts, m_pop)
            self.motif_states.append(motif_state)
        self.motif_states.sort(key = lambda x : x.population, reverse=True)




