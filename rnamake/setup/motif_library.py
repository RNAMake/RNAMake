import os
import rnamake.settings as settings
import rnamake.motif_factory as motif_factory
import rnamake.motif_type as motif_type
import rnamake.util as util

class MotifLibrary(object):
    def __init__(self, libtype=None, libdir=None, libfile=None):
        self.motif_paths = {}
        self.mtype = motif_type.UNKNOWN
        self.motif_dict = {}
        self._parse_args(libtype, libdir, libfile)

    def get_motif(self, mname):
        if mname not in self.motif_paths:
            raise ValueError("unknown motif: " + mname + "cannot get it")

        if mname not in self.motif_dict:
            self.motif_dict[mname] = motif_factory.factory.motif_from_file(self.motif_paths[mname])
            if self.motif_dict[mname] is None:
                return None
            self.motif_dict[mname].mtype = self.mtype

        return self.motif_dict[mname].copy()

    def load_all(self, limit=999999):
        for i, mname in enumerate(self.motif_paths):
            self.get_motif(mname)
            if i > limit:
                return

    def motifs(self):
        return self.motif_dict.values()

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


def unique_twoway_lib():
    path = settings.MOTIF_DIRS + "two_ways/unique_7.dat"
    mlib = MotifLibrary(libfile=path)
    mlib.load_all()
    return mlib

def ideal_helix_lib():
    mlib = MotifLibrary(motif_type.HELIX)
    mlib.get_motif("HELIX.IDEAL")
    for i in range(1,21):
        mlib.get_motif("HELIX.IDEAL."+str(i))
    return mlib


lib_paths = {
    motif_type.TWOWAY   : "two_ways",
    motif_type.HAIRPIN  : "hairpins",
    motif_type.HELIX    : "helices",
    motif_type.NWAY     : "junctions",
    motif_type.TCONTACT : "tertiary_contacts"
}
