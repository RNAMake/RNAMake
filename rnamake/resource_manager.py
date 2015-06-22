import motif
import motif_factory
import sqlite_library
import motif_type
import settings
import util

class ResourceManager(object):
    def __init__(self):
        self.mlibs = {}
        self.extra_motifs = {}

        self.ms_libs = {}
        self.extra_ms = {}

        for k in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            self.mlibs[k] = sqlite_library.MotifSqliteLibrary(k)

        for k in sqlite_library.MotifStateSqliteLibrary.get_libnames().keys():
            self.ms_libs[k] = sqlite_library.MotifStateSqliteLibrary(k)

    def get_motif(self, mname, end_index=None, end_name=None):
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

    def add_motif(self, path):
        name = util.filename(path)
        self.extra_motifs[name] = motif.Motif(path)


manager = ResourceManager()

