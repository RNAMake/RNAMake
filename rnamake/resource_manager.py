import motif
import motif_library
import motif_library_sqlite
import motif_type
import motif_tree_state
import motif_tree_state_library
import motif_scorer
import settings
import util

class ResourceManager(object):
    def __init__(self):
        self.mlibs = {}
        self.extra_motifs = {}

        self.mts_libs = {}
        self.extra_mts = {}

        for libname in motif_library_sqlite.libnames.keys():
            # catch unimplemented libraries
            mlib = motif_library_sqlite.MotifLibrarySqlite(libname=libname)
            self.mlibs[libname] = mlib

        path = settings.RESOURCES_PATH + "motif_libraries/bp_steps.db"
        mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)
        mlib.mtype = motif_type.HELIX

        path = settings.RESOURCES_PATH + "motif_libraries/seq_and_ss.db"
        mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)
        mlib.mtype = motif_type.HELIX

        self.mlibs['UNKNOWN'] = mlib

        self._setup_mts_libs()
        #self.mt = motif_tree.MotifTree()
        #self.mt2 = motif_tree.MotifTree()
        #self.mt.add_motif(self.get_motif("HELIX.IDEAL.6"), end_index=1, end_flip=0)
        #self.mt.level +=1

    def _setup_mts_libs(self):
        path = settings.UNITTEST_PATH
        self.mts_libs['TWOWAY']   = motif_tree_state_library.MotifTreeStateLibrary(libpath=path+"test_twoway.new.me",new=1)
        self.mts_libs['HELIX' ]   = motif_tree_state_library.MotifTreeStateLibrary(libpath=path+"test_helix.new.me",new=1)
        self.mts_libs['BP_STEPS'] = motif_tree_state_library.MotifTreeStateLibrary(libpath=settings.RESOURCES_PATH+"prediction/all.new.me",new=1)

    def get_motif(self, mname, end_index=None, end_name=None):
        for mlib in self.mlibs.itervalues():
            if mname in mlib:
                return mlib.get_motif(mname)

        """if mname in self.extra_motifs:
            if end_index is None and end_name is None:
                return self.extra_motifs[mname]

            m = self.extra_motifs[mname]
            return self._prep_extra_motif_for_asssembly(m, end_index, end_name)"""


        raise ValueError("cannot find " + mname)

    def get_state(self, mts_name):
        for mts_lib in self.mts_libs.itervalues():
            mts = mts_lib.get_state(mts_name)
            if mts is not None:
                return mts

        if mts_name in self.extra_mts:
            return self.extra_mts[mts_name]

        raise ValueError("cannot find mts: "+ mts_name)

    def add_lib_path(self, path):
        self.mlibs[path] = motif_library.MotifLibrary(libdir=path)

    def add_motif(self, path):
        name = util.filename(path)
        self.extra_motifs[name] = motif.Motif(path)

    def _prep_extra_motif_for_asssembly(self, m, end_index=None, end_name=None):
        if end_name is None and end_index is None:
            raise ValueError("cannot call _prep_extra_motif_for_asssembly without end_name" +\
                             "or end_index" )

        if end_name is not None:
            end = m.get_basepair_by_name(end_name)
            end_index = m.ends.index(end)

        node = self.mt.add_motif(m, end_index=end_index)
        if node.flip:
            node.motif.ends[end_index].flipped=0

        if node is None:
            raise ValueError("cannot prep motif for assembly motif:" + m.name + \
                             " end_index: " + str(end_index))

        h_m = self.get_motif("HELIX.IDEAL.3")
        for e in node.available_ends():
            n = self.mt.add_motif(h_m, end_index=1, end_flip=0, parent=node)
            if n is None:
                e.flip()
            else:
                self.mt.remove_node(self.mt.last_node)

        m_copy = self.mt.nodes[2].motif
        self.mt.remove_node_level()
        self._build_extra_mts(m_copy, end_index)

        return m_copy

    def _build_extra_mts(self, m, end_index):
        self.mt2.add_motif(m, end_index=end_index, end_flip=0)
        m_copy = self.mt2.nodes[1].motif
        self.mt2.remove_node_level()

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

        name = m_copy.name + "-" + str(end_index)
        mts = motif_tree_state.MotifTreeState(name, end_index,
                                              len(m.residues()), score, beads,
                                              ends, 0, m_copy.to_str())

        self.extra_mts[name] = mts



manager = ResourceManager()

