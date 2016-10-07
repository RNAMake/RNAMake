import motif
import motif_factory
import sqlite_library
import motif_type
import motif_ensemble
import settings
import util
import exceptions

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
        motifs[0].new_res_uuids()
        return motifs[0]

    def get_multi(self, **options):
        motifs = self._find_motifs(**options)
        for m in motifs:
            m.new_res_uuids()
        return motifs

class ResourceManager(object):
    def __init__(self):
        self.mlibs,    self.added_motifs = {}, MotifLibrary()
        self.ms_libs,  self.extra_ms = {}, {}
        self.me_libs,  self.extra_me = {}, {}
        self.mse_libs, self.extra_mse = {}, {}
        exclude = ['all_bp_steps']

        for k in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            if k in exclude:
                continue
            self.mlibs[k] = sqlite_library.MotifSqliteLibrary(k)

        for k in sqlite_library.MotifStateSqliteLibrary.get_libnames().keys():
            if k in exclude:
                continue
            self.ms_libs[k] = sqlite_library.MotifStateSqliteLibrary(k)

        for k in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
            if k in exclude:
                continue
            self.me_libs[k] = sqlite_library.MotifEnsembleSqliteLibrary(k)

        for k in sqlite_library.MotifStateEnsembleSqliteLibrary.get_libnames().keys():
            if k in exclude:
                continue
            self.mse_libs[k] = sqlite_library.MotifStateEnsembleSqliteLibrary(k)


        #hack for now
        self.add_motif(settings.MOTIF_DIRS + "/extras/GAAA_tetraloop")
        self.add_motif(settings.MOTIF_DIRS + "/extras/GGAA_tetraloop")

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

        raise exceptions.ResourceManagerException(
            "cannot find motif: " + self._args_to_str(options))

    def contains_motif(self, **options):
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

        raise exceptions.ResourceManagerException(
            "cannot find motif state: "+ self._args_to_str(options))

    def get_motif_ensemble(self, ss_id):
        for me_lib in self.me_libs.itervalues():
            if me_lib.contains(name=ss_id):
                return me_lib.get(name=ss_id)

        raise ValueError("could not find motif ensemble with id: "+ss_id)

    def get_motif_state_ensemble(self, **options):
        for mse_lib in self.mse_libs.itervalues():
            if mse_lib.contains(**options):
                return mse_lib.get(**options)

        raise ValueError("could not find motif state ensemble with options:"+\
                         self._args_to_str(options))

    def add_motif(self, path=None, motif=None, name=None, include_protein=0,
                  align=1):
        if path:
            m = motif_factory.factory.motif_from_file(path,
                                                      include_protein=include_protein)
        else:
            m = motif

        if name is not None:
            m.name = name

        if len(m.ends) == 0:
            self.added_motifs.add_motif(m)
            return

        motifs = []
        end_ids = {}
        for i in range(len(m.ends)):
            m_added = motif_factory.factory.can_align_motif_to_end(m, i)
            if m_added is None:
                continue
            if not align:
                motif_factory.factory.standardize_motif(m_added)
                motifs.append(m)
                end_ids[m_added.ends[0].uuid] = m_added.end_ids[0]
                break

            m_added = motif_factory.factory.align_motif_to_common_frame(m_added, i)
            if m_added is None:
                continue
            motifs.append(m_added)
            end_ids[m_added.ends[0].uuid] = m_added.end_ids[0]

        #fix minor changes between ids
        for m in motifs:
            for i, end in enumerate(m.ends):
                try :
                    end_id = end_ids[end.uuid]
                    m.end_ids[i] = end_id
                except:
                    pass
            self.added_motifs.add_motif(m)

    def register_motif(self, m):
        if m.name == "":
            raise exceptions.ResourceManagerException(
                "attempted to register motif with no name this will make it "
                "extremely unlikely you will be able to retrieve it properly!")

        self.added_motifs.add_motif(m)

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


manager = ResourceManager()

