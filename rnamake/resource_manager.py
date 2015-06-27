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

        for k in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            self.mlibs[k] = sqlite_library.MotifSqliteLibrary(k)

        for k in sqlite_library.MotifStateSqliteLibrary.get_libnames().keys():
            self.ms_libs[k] = sqlite_library.MotifStateSqliteLibrary(k)

        for k in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
            self.me_libs[k] = sqlite_library.MotifEnsembleSqliteLibrary(k)

        for k in sqlite_library.MotifStateEnsembleSqliteLibrary.get_libnames().keys():
            self.mse_libs[k] = sqlite_library.MotifStateEnsembleSqliteLibrary(k)

    def get_motif(self, mname=None, ss_id=None):

        if mname is not None:
            for mlib in self.mlibs.itervalues():
                if mlib.contains(mname):
                    return mlib.get(mname)

        raise ValueError("cannot find " + mname)

    def get_state(self, ms_name):
        for ms_lib in self.ms_libs.itervalues():
            if ms_lib.contains(ms_name):
                return ms_lib.get(ms_name)

        if ms_name in self.extra_ms:
            return self.extra_ms[ms_name]

        raise ValueError("cannot find mts: "+ mts_name)

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
        name = util.filename(path)
        self.extra_motifs[name] = motif.Motif(path)


manager = ResourceManager()

