import motif
import motif_factory
import sqlite_library
import motif_type
import motif_ensemble
import settings
import util
import exceptions
import residue_type

class MotifLibrary(object):
    def __init__(self):
        self.motifs = []

    def add_motif(self, m):
        self.motifs.append(m)

    def _find_motifs(self, **options):
        motifs = []
        for m in self.motifs:
            if 'name' in options:
                if options['name'] != m.name:
                    continue
            if 'end_name' in options:
                if options['end_name'] != m.get_end(0).name:
                    continue
            if 'end_id' in options:
                if options['end_id'] != m.get_end_id(0):
                    continue
            if m not in motifs:
                motifs.append(m)
        return motifs

    def contains(self, **options):
        motifs = self._find_motifs(**options)
        if len(motifs) > 0:
            return 1
        else:
            return 0

    def get(self, **options):
        motifs = self._find_motifs(**options)
        return motif.Motif.copy(motifs[0], new_uuid=1)

    def get_multi(self, **options):
        motifs = self._find_motifs(**options)
        return [motif.Motif.copy(m, new_uuid=1) for m in motifs]


class ResourceManager(object):

    __slots__ = [
        "_mlibs",
        "_me_libs",
        "_added_motifs",
        "_mf",
        "_rts"
    ]

    def __init__(self):
        self._mlibs = {}
        self._me_libs = {}
        self._added_motifs = MotifLibrary()
        self._rts = residue_type.ResidueTypeSet()
        self._mf = motif_factory.MotifFactory(self._rts)

        for k in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            self._mlibs[k] = sqlite_library.MotifSqliteLibrary(k, self._rts)

        for k in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
            self._me_libs[k] = sqlite_library.MotifEnsembleSqliteLibrary(k, self._rts)

        # hack for now
        self.add_motif_from_file(settings.MOTIF_DIRS + "/extras/GAAA_tetraloop")
        self.add_motif_from_file(settings.MOTIF_DIRS + "/extras/GGAA_tetraloop")

    def get_motif(self, **options):
        for mlib in self._mlibs.itervalues():
            if mlib.contains(**options):
                return mlib.get(**options)

        if self._added_motifs.contains(**options):
            return self._added_motifs.get(**options)

        raise exceptions.ResourceManagerException(
            "cannot find motif: " + self._args_to_str(options))

    def contains_motif(self, **options):
        for mlib in self._mlibs.itervalues():
            if mlib.contains(**options):
                return 1

        if self._added_motifs.contains(**options):
            return 1

        return 0

    def get_bp_step(self, end_id):
        m = self._mlibs['new_bp_steps'].get(end_id=end_id)
        return m

    def get_motif_multi(self, **options):
        for mlib in self.mlibs.itervalues():
            if mlib.contains(**options):
                return mlib.get_multi(**options)

        if self.added_motifs.contains(**options):
            return self.added_motifs.get_multi(**options)

        raise ValueError("cannot find motif")

    def get_state(self, **options):
        for mlib in self._mlibs.itervalues():
            if mlib.contains(**options):
                return mlib.get(**options).get_state()

        if self._added_motifs.contains(**options):
            return self._added_motifs.get(**options).get_state()

        raise exceptions.ResourceManagerException(
            "cannot find motif state: "+ self._args_to_str(options))

    def get_motif_ensemble(self, ss_id):
        new_me = None
        for me_lib in self._me_libs.itervalues():
            if me_lib.contains(name=ss_id):
                new_me = me_lib.get(name=ss_id)
                break

        if new_me is not None:
            for mem in new_me:
                self.add_motif(mem.motif, mem.motif.mtype, mem.motif.name)
            return new_me

        raise ValueError("could not find motif ensemble with id: "+ss_id)

    def get_motif_state_ensemble(self, **options):
        for mse_lib in self.mse_libs.itervalues():
            if mse_lib.contains(**options):
                return mse_lib.get(**options)

        raise ValueError("could not find motif state ensemble with options:"+\
                         self._args_to_str(options))

    def add_motif_from_file(self, path=None, name=None, include_protein=0,
                  align=1):

        motifs = self._mf.motifs_from_file(path, include_protein=include_protein)

        for m in motifs:
            self._added_motifs.add_motif(m)

    def add_motif(self, m, mtype=motif_type.UNKNOWN, mname="unknown"):
        motifs = self._mf.motifs_from_rstruc(m, mtype, mname)

        for m in motifs:
            self._added_motifs.add_motif(m)

    def register_motif(self, m):
        if m.name == "":
            raise exceptions.ResourceManagerException(
                "attempted to register motif with no name this will make it "
                "extremely unlikely you will be able to retrieve it properly!")

        self.added_motifs.add_motif(m.copy())

    def register_extra_motif_ensembles(self, f_name):
        f = open(f_name)
        lines = f.readlines()
        f.close()

        for l in lines:
            spl = l.split("!!")
            self.extra_me[spl[0]] = motif_ensemble.str_to_motif_ensemble(spl[1])

    def has_supplied_motif_ensemble(self, m_name, end_name):
        key = m_name + "-" + end_name

        if key in self.extra_me:
            return 1
        else:
            return 0

    def get_supplied_motif_ensemble(self, m_name, end_name):
        key = m_name + "-" + end_name

        if key not in self.extra_me:
            raise ValueError("no supplied motif ensemble exists")

        return self.extra_me[key]

    def _args_to_str(self, options):
        s = ""
        for k, v in options.iteritems():
            s += k + " = " + v + ","
        return s

    def ms_to_motif(self, ms):
        m = self.get_motif(name=ms.name, end_name=ms.get_end(0).name)
        pass

    def get_motif_with_new_alignment(self, m, ei):
        m_copy = motif.Motif.copy(m)
        m_copy.get_end(ei).flip()

        new_m = self.get_motif(name=m_copy.name, end_name=m.get_end(ei).name)
        m_aligned = motif.get_aligned_motif(m_copy.get_end(ei), new_m.get_end(0), new_m)
        return m_aligned

    @property
    def residue_type_set(self):
        return self._rts



