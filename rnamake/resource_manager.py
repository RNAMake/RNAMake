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

        #if 'name' in options:
        #    motifs = self.get_motif_multi(name=options['name'])
        #    s = "cannot find motif: " + self._args_to_str(options)  + "\n"
        #    s += "here are the available options: \n"
        #    for m in motifs:
        #        s += "name: %s, end_id: %s, end_name: %s\n" % (m.name, m.end_ids[0],
        #                                                       m.ends[0].name())
        #    raise ValueError(s)

        raise ValueError("cannot find motif: " + self._args_to_str(options))

    def motif_exists(self, **options):
        for mlib in self.mlibs.itervalues():
            if mlib.contains(**options):
                return 1

        if self.added_motifs.contains(**options):
            return 1

        return 0

    def _args_to_str(self, options):
        s = ""
        for k, v in options.iteritems():
            s += k + " = " + v + ","
        return s

    def get_motif_multi(self, **options):
        for mlib in self.mlibs.itervalues():
            if mlib.contains(**options):
                return mlib.get_multi(**options)

        if self.added_motifs.contains(**options):
            return self.added_motifs.get_multi(**options)

        raise ValueError("cannot find motif")

    def get_state(self, **options):
        for me_lib in self.ms_libs.values():
            if me_lib.contains(**options):
                return me_lib.get(**options)

        if self.added_motifs.contains(**options):
            return self.added_motifs.get(**options).get_state()

        raise ValueError("cannot find mts: "+ ms_name)

    def get_motif_ensemble(self, ss_id):
        for me_lib in self.me_libs.itervalues():
            if me_lib.contains(ss_id):
                return me_lib.get(ss_id)

        raise ValueError("could not find motif ensemble with id: "+ss_id)

    def get_motif_state_ensemble(self, **options):
        for mse_lib in self.mse_libs.itervalues():
            if mse_lib.contains(**options):
                return mse_lib.get(**options)

        raise ValueError("could not find motif state ensemble with options:"+\
                         self._args_to_str(options))

    def add_motif(self, path=None, motif=None):
        if path:
            m = motif_factory.factory.motif_from_file(path)
        else:
            m = motif
        motifs = []
        end_ids = {}
        for i in range(len(m.ends)):
            m_added = motif_factory.factory.can_align_motif_to_end(m, i)
            if m_added is None:
                continue
            m_added = motif_factory.factory.align_motif_to_common_frame(m_added, i)
            if m_added is None:
                continue
            motifs.append(m_added)
            end_ids[m_added.ends[0].uuid] = m_added.end_ids[0]

        #fix minor changes between ids
        for m in motifs:
            for i, end in enumerate(m.ends):
                end_id = end_ids[end.uuid]
                m.end_ids[i] = end_id
            self.added_motifs.add_motif(m)

    def register_motif(self, m):
        self.added_motifs.add_motif(m)

manager = ResourceManager()

