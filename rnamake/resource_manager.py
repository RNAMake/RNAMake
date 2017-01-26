import motif
import motif_factory
import sqlite_library
import motif_type
import motif_ensemble
import settings
import util
import exceptions
import residue_type
import residue
import chain
import structure
import basepair

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
        "_ms_libs",
        "_added_motifs",
        "_mf",
        "_rts"
    ]

    def __init__(self):
        self._mlibs = {}
        self._me_libs = {}
        self._ms_libs = {}
        self._added_motifs = MotifLibrary()
        self._rts = residue_type.ResidueTypeSet()
        self._mf = motif_factory.MotifFactory(self._rts)

        for k in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            self._mlibs[k] = sqlite_library.MotifSqliteLibrary(k, self._rts)

        for k in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
            self._me_libs[k] = sqlite_library.MotifEnsembleSqliteLibrary(k, self._rts)

        for k in sqlite_library.MotifStateSqliteLibrary.get_libnames().keys():
            self._ms_libs[k] = sqlite_library.MotifStateSqliteLibrary(k)

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

    def get_bp_step(self, end_id, org_m=None):
        m = self._mlibs['bp_steps'].get(end_id=end_id)
        if org_m is None:
            return m

        new_chains = []
        for i, c1 in enumerate(m.iter_chains()):
            c2 = org_m.get_chain(i)
            new_res = []
            for j, r1 in enumerate(c1):
                new_r = residue.Residue.copy(r1, given_uuid=c2.get_residue(j).uuid)
                new_res.append(new_r)
            c = chain.Chain(new_res)
            new_chains.append(c)
        s = structure.Structure(new_chains)
        new_bps = []
        for i, bp in enumerate(m.iter_basepairs()):
            bp_res = m.get_bp_res(bp)
            r_pos_1 = m.get_res_index(bp_res[0])
            r_pos_2 = m.get_res_index(bp_res[1])
            other_res1 = org_m.get_residue(index=r_pos_1)
            other_res2 = org_m.get_residue(index=r_pos_2)
            bp2 = org_m.get_basepair(uuid1=other_res1.uuid, uuid2=other_res2.uuid)
            new_bp = basepair.Basepair.copy_with_new_uuids(
                            bp, bp2.res1_uuid, bp2.res2_uuid, bp2.uuid)
            new_bps.append(new_bp)

        new_ends = []
        for end in m.iter_ends():
            for new_bp in new_bps:
                if end.name == new_bp.name:
                    new_ends.append(new_bp)
                    break

        end_ids = [ m.get_end_id(i) for i in range(m.num_ends()) ]
        new_m = motif.Motif(s, new_bps, new_ends, end_ids, m.name, m.mtype,
                            m.score, m.block_end_add, m.dot_bracket, [], org_m.uuid)
        return new_m

    def get_motif_multi(self, **options):
        for mlib in self.mlibs.itervalues():
            if mlib.contains(**options):
                return mlib.get_multi(**options)

        if self.added_motifs.contains(**options):
            return self.added_motifs.get_multi(**options)

        raise ValueError("cannot find motif")

    def get_state(self, **options):

        for ms_lib in self._ms_libs.itervalues():
            if ms_lib.contains(**options):
                return ms_lib.get(**options)

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

    def add_motif_from_file(self, path=None, name=None, include_protein=0):

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

        self._added_motifs.add_motif(m.copy())

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

    @property
    def motif_factory(self):
        return self._mf

    def get_mlib(self, name):
        if name in self._mlibs:
            return self._mlibs[name]
        else:
            raise exceptions.ResourceManagerException("cannot find motif lib: " + name)

    def get_ms_lib(self, name):
        if name in self._ms_libs:
            return self._ms_libs[name]
        else:
            raise exceptions.ResourceManagerException("cannot find ms lib: " + name)