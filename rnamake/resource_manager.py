from . import motif_library
from . import motif_type

class ResourceManager(object):
    def __init__(self):
        self.mlibs = {}
        self.mts_libs = {}

        for mtype in motif_library.lib_paths.iterkeys():
            #catch uninmplemented libraries
            try:
                mlib = motif_library.MotifLibrary(mtype)
                self.mlibs[motif_type.type_to_str(mtype)] = mlib
            except:
                pass

    def get_motif(self, mname):
        for mlib in self.mlibs.itervalues():
            if mname in mlib:
                return mlib.get_motif(mname)

        raise ValueError("cannot find " + mname)


