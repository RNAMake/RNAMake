import motif
import motif_factory
import sqlite_library
import motif_type
import settings
import util

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
                if options['end_name'] != m.ends[0].name():
                    continue
            if 'end_id' in options:
                if options['end_id'] != m.end_ids[0]:
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
        return motifs[0]

    def get_multi(self, **options):
        return self._find_motifs(**options)

class ResourceManager(object):
    def __init__(self):
        self.mlibs,    self.added_motifs = {}, MotifLibrary()
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

    def get_motif(self, **options):
        for mlib in self.mlibs.itervalues():
            if mlib.contains(**options):
                return mlib.get(**options)

        if self.added_motifs.contains(**options):
            return self.added_motifs.get(**options)

        raise ValueError("cannot find motif")

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
            self.added_motifs.add_motif(m_added)


manager = ResourceManager()

