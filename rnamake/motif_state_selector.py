import motif_type
import settings
import resource_manager
import graph
import exceptions
import motif_state
import util


class MotifStateSelector(object):
    __slots__ = [
        '_rm',
        '_graph'
    ]

    def __init__(self, rm):
        self._rm = resource_manager.ResourceManager()
        self._graph = graph.GraphDynamic()

    def __get_motif_states_by_name(self, name):
        try:
            mlib = self._rm.get_ms_lib(name)
            mlib.load_all()
            return motif_state.MotifLibrary(name, mlib.all())
        except exceptions.ResourceManagerException:
            pass

        if name == 'unique_twoway':
            path = settings.RESOURCES_PATH+"/motif_lists/unique_twoways.mlist"
            lines = util.get_file_contents(path)
            motif_states = []
            for l in lines:
                spl = l.split("|")
                m_name, end_name = spl[0].split(",")
                try:
                    m_state = self._rm.get_state(name=m_name, end_name=end_name)
                except:
                    continue
                motif_states.append(m_state)
            return motif_state.MotifLibrary(name, motif_states)

        if name == 'ideal_helices_min':
            mlib = self._rm.get_ms_lib("ideal_helices")
            mlib.load_all()
            motif_states = []
            for m in mlib.all():
                if m.name == "HELIX.IDEAL":
                    continue
                motif_states.append(m)
            return motif_state.MotifLibrary(name, motif_states)

        raise ValueError("motif library: " + name  + " is not recognized")

    def add(self, name=None, mlib=None, max_uses=1000, required_uses=0):
        if name is not None:
            named_mlib = self.__get_motif_states_by_name(name)
            d = MotifLibraryData(named_mlib, max_uses, required_uses)
            self._graph.add_data(d, orphan=1)

        if mlib is not None:
            d = MotifLibraryData(mlib, max_uses, required_uses)
            self._graph.add_data(d, orphan=1)

    def connect(self, name_i, name_j):
        i, j = -1, -1
        for n in self._graph:
            if n.data.mlib.name == name_i:
                i = n.index
            if n.data.mlib.name == name_j:
                j = n.index

        if i == -1 or j == -1:
            raise ValueError("could not find a node with that name")

        self._graph.connect(i, j)

    def get_connected_libs(self, ntype):
        connected_indexes = []
        if ntype != -1:
            connections = self._graph.get_node(ntype).connections
            for c in connections:
                connected_indexes.append(c.partner(ntype).index)
        else:
            connected_indexes = [0]
        return connected_indexes

    def get_motif_states(self, ntype):
        return self._graph.get_node(ntype).data.mlib.motifs


class MotifLibraryData(object):
    def __init__(self, mlib, max_uses, required_uses):
        self.mlib = mlib


def default_selector(rm):
    selector = MotifStateSelector(rm)
    selector.add('ideal_helices_min')
    selector.add('unique_twoway')
    selector.connect('ideal_helices_min', 'unique_twoway')
    return selector










