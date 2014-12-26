import os
from . import settings
from . import motif_type
from . import motif
from . import util

class MotifLibrary(object):
    def __init__(self, libtype=None, libdir=None, libfile=None):
        self.motif_paths = {}
        self.mtype = motif_type.UNKNOWN
        self.motifs = {}
        self._parse_args(libtype, libdir, libfile)

    def get_motif(self, mname):
        if mname not in self.motif_paths:
            raise ValueError("unknown motif: " + mname + "cannot get it")

        if mname not in self.motifs:
            self.motifs[mname] = motif.Motif(self.motif_paths[mname])
            self.motifs[mname].mtype = self.mtype

        return self.motifs[mname].copy()

    def _parse_args(self, libtype, libdir, libfile):
        if libtype is None and libdir is None and libfile is None:
            raise ValueError("please specify either the type, dir or lib " +\
                             "cannot initiate without one")

        if   libtype is not None:
            if libtype not in lib_paths:
                raise ValueError("libtype " + libtype + "is not recognized")
            self.mtype = libtype
            type_dir = settings.MOTIF_DIRS + lib_paths[libtype]
            self._get_motifs_from_dir(type_dir)
        elif libdir is not None:
            self._get_motifs_from_dir(libdir)
        else:
            self._get_motifs_from_file(libfile)

    def _get_motifs_from_dir(self, libdir):

        for f in os.listdir(libdir):
            mdir = libdir + "/" + f
            if not os.path.isdir(mdir) or f[:1] == ".":
                continue
            if not os.path.isfile(mdir + "/ref_frames.dat"):
                continue
            self.motif_paths[f] = mdir

        if len(self.motif_paths.keys()) == 0:
            raise ValueError("no motifs were loaded from " + libdir)

    def _get_motifs_from_file(self, libfile):
        try:
            f = open(libfile)
            lines = f.readlines()
            f.close()
        except IOError:
            raise IOError("motif list file: " + libfile + " does not exist \
                          please check the path")

        libdir = util.base_dir(libfile)
        for l in lines:
            spl = l.split()
            self.motif_paths[spl[0]] = libdir + "/" + spl[0]

    def __contains__(self, mname):
        return mname in self.motif_paths


lib_paths = {
    motif_type.TWOWAY   : "two_ways",
    motif_type.HAIRPIN  : "hairpins",
    motif_type.HELIX    : "helices",
    motif_type.NWAY     : "junctions"
}
