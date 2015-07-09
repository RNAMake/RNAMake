import motif
import motif_factory
import sqlite_library
import motif_type
import settings
import util

class ResourceManager(object):
    def __init__(self):
        self.mlibs,    self.extra_motifs = {}, {}
        self.ms_libs,  self.extra_ms = {}, {}
        self.me_libs,  self.extra_me = {}, {}
        self.mse_libs, self.extra_mse = {}, {}
        self.ss_libs,  self.extra_ss = {}, {}

        for k in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            self.mlibs[k] = sqlite_library.MotifSqliteLibrary(k)

        for k in sqlite_library.MotifStateSqliteLibrary.get_libnames().keys():
            self.ms_libs[k] = sqlite_library.MotifStateSqliteLibrary(k)

        for k in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
            self.me_libs[k] = sqlite_library.MotifEnsembleSqliteLibrary(k)

        for k in sqlite_library.MotifStateEnsembleSqliteLibrary.get_libnames().keys():
            self.mse_libs[k] = sqlite_library.MotifStateEnsembleSqliteLibrary(k)

        for k in sqlite_library.MotifSSIDSqliteLibrary.get_libnames().keys():
            self.ss_libs[k] = sqlite_library.MotifSSIDSqliteLibrary(k)

    def get_motif(self, mname, end_name=None):
        if end_name is not None:
            mname = mname + "-" + end_name

        for mlib in self.mlibs.itervalues():
            if mlib.contains(mname):
                return mlib.get(mname)

        for mlib in self.ss_libs.itervalues():
             if mlib.contains(mname):
                 return mlib.get(mname)

        if mname in self.extra_motifs:
            return self.extra_motifs[mname]


        raise ValueError("cannot find " + mname)

    def get_state(self, ms_name, end_name=None):
        if end_name is not None:
            ms_name = ms_name + "-" + end_name

        for ms_lib in self.ms_libs.itervalues():
            if ms_lib.contains(ms_name):
                return ms_lib.get(ms_name)

        if ms_name in self.extra_motifs:
            return self.extra_motifs[ms_name].get_state()

        raise ValueError("cannot find mts: "+ ms_name)

    def get_motif_ensemble(self, ss_id):
        for me_lib in self.me_libs.itervalues():
            if me_lib.contains(ss_id):
                return me_lib.get(ss_id)

        raise ValueError("could not find motif ensemble with id: "+ss_id)

    def get_motif_state_ensemble(self, ss_id):
        for mse_lib in self.mse_libs.itervalues():
            if mse_lib.contains(ss_id):
                return mse_lib.get(ss_id)

        raise ValueError("could not find motif state ensemble with id: "+ss_id)

    def add_motif(self, path):
        m = motif_factory.factory.motif_from_file(path)
        for i in range(len(m.ends)):
            m_added = motif_factory.factory.can_align_motif_to_end(m, i)
            if m_added is None:
                continue
            m_added = motif_factory.factory.align_motif_to_common_frame(m_added, i)
            if m_added is None:
                continue
            m_added.name = m_added.name + "-" + m_added.ends[0].name()
            self.extra_motifs[m_added.name] = m_added


manager = ResourceManager()

