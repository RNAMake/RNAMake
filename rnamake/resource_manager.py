import motif_library
import motif_library_sqlite
import motif_type
import motif_tree_state

class ResourceManager(object):
    def __init__(self):
        self.mlibs = {}
        self.mts_libs = {}

        for mtype in motif_library.lib_paths.iterkeys():
            #catch uninmplemented libraries
            try:
                mlib = motif_library_sqlite.MotifLibrarySqlite(mtype)
                self.mlibs[motif_type.type_to_str(mtype)] = mlib
            except:
                pass

        for mtype in motif_library.lib_paths.iterkeys():
            #catch unimplemented mts libraries
            try:
                mts_lib = motif_tree_state.MotifTreeStateLibrary(mtype)
                self.mts_libs[motif_type.type_to_str(mtype)] = mts_lib
            except:
                pass

    def get_motif(self, mname):
        for mlib in self.mlibs.itervalues():
            if mname in mlib:
                return mlib.get_motif(mname)

        raise ValueError("cannot find " + mname)

    def add_lib_path(self, path):
        self.mlibs[path] = motif_library.MotifLibrary(libdir=path)



