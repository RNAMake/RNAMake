import motif_library
import rnamake.motif_factory as motif_factory
import rnamake.motif_ensemble as motif_ensemble
import rnamake.motif as motif
import rnamake.motif_type as motif_type
import rnamake.util as util
import rnamake.cluster as cluster
import rnamake.settings as settings
import rnamake.sqlite_library as sqlite_library
import math
import subprocess

def test_align():
    mlib = sqlite_library.MotifLibrarySqlite(libname="ideal_helices")
    m1 = mlib.get_motif("HELIX.IDEAL.12")
    m2 = m1.copy()

    aligned = motif.get_aligned_motif(m1.ends[1], m2.ends[0], m2)

    m1.to_pdb("test_1.pdb")
    aligned.to_pdb("test_2.pdb")

def setup_start_motif():
    start_path = settings.RESOURCES_PATH + "start"
    m = motif_factory.factory.motif_from_file(start_path)
    s = m.to_str()
    m_path = settings.MOTIF_DIRS + "ref.motif"
    print m_path
    f = open(m_path, "w")
    f.write(s)
    f.close()

    m2 = motif.file_to_motif(m_path)

class SSandSeqCluster(object):

    def __init__(self, end_id):
        self.end_id = end_id
        self.motif_and_ends = []

    def motif_matches_end(self, m, ei):
        class MotifandEnd(object):
            def __init__(self, m, ei):
                self.motif, self.end_index = m, ei

        if self.end_id == m.end_ids[ei]:
            self.motif_and_ends.append(MotifandEnd(m, ei))
            return 1
        else:
            return 0

    def motif_matches(self, m):
        for i, end_id in enumerate(m.end_ids):
            r = self.motif_matches_end(m, i)
            if r == 1:
                return 1
        return 0


class BuildSqliteLibraries(object):

    def build_ideal_helices(self):
        mlib = motif_library.ideal_helix_lib()
        aligned_motifs = []
        names = []
        path = settings.RESOURCES_PATH +"/motif_libraries_new/ideal_helices.db"
        data = []
        keys = ['data', 'name', 'end_name', 'end_id', 'id']
        for i, m in enumerate(mlib.motifs()):
            aligned_motif = motif_factory.factory.align_motif_to_common_frame(m, 1)
            aligned_motifs.append(aligned_motif)
            names.append(aligned_motif.name)
            data.append([aligned_motif.to_str(), aligned_motif.name,
                         aligned_motif.ends[0].name(), aligned_motif.end_ids[0], i])

        sqlite_library.build_sqlite_library_2(path, data, keys, 'id')
        mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        m = mlib.get(name="HELIX.IDEAL.6")
        s = m.to_str()
        path = settings.RESOURCES_PATH + "/motifs/base.motif"
        f = open(path, "w")
        f.write(s)
        f.close()

    def build_basic_libraries(self):

        types = [motif_type.TWOWAY, motif_type.NWAY, motif_type.HAIRPIN,
                 motif_type.TCONTACT]

        for t in types:
            count = 0
            mlib = motif_library.MotifLibrary(t)
            mlib.load_all()

            succeses = []
            for k, m in enumerate(mlib.motifs()):
                if m is None:
                    continue

                if t != motif_type.HAIRPIN and len(m.ends) == 1:
                    continue
                elif len(m.ends) == 0:
                    continue

                for ei in range(len(m.ends)):
                    m_added = motif_factory.factory.can_align_motif_to_end(m, ei)
                    if m_added is None:
                        continue
                    else:
                        succeses.append([m_added, ei])

            data = []
            keys = ['data', 'name', 'end_name', 'end_id', 'id']
            for i, s in enumerate(succeses):
                m, ei = s
                m_added = motif_factory.factory.align_motif_to_common_frame(m, ei)
                data.append([m_added.to_str(), m_added.name,
                            m_added.ends[0].name(), m_added.end_ids[0], i])

            path = settings.RESOURCES_PATH +"/motif_libraries_new/"+\
                   motif_type.type_to_str(t).lower()+".db"
            sqlite_library.build_sqlite_library_2(path, data, keys, 'id')

    def build_helix_ensembles(self):
        helix_mlib = motif_library.MotifLibrary(motif_type.HELIX)
        helix_mlib.load_all()

        clusters = []
        steps = []

        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

        for m in helix_mlib.motifs():
            spl = m.name.split(".")
            if spl[1] == "IDEAL" or spl[1] == "LE":
                continue
            for i in range(len(m.basepairs)-1):
                if m.basepairs[i].bp_type != "cW-W" or m.basepairs[i+1].bp_type != "cW-W":
                    continue
                bps = [m.basepairs[i], m.basepairs[i+1]]
                res = []
                for bp in bps:
                    for r in bp.residues():
                        if r not in res:
                            res.append(r)
                if len(res) != 4:
                    continue
                m_bps = motif_factory.factory.motif_from_bps(bps)
                if len(m_bps.end_ids[0]) != 11:
                    continue
                matched = 0
                for c in clusters:
                    if c.motif_matches(m_bps):
                        matched = 1

                if not matched:
                    clusters.append(SSandSeqCluster(m_bps.end_ids[0]))
                    clusters[-1].motif_matches(m_bps)
                    if m_bps.end_ids[0] != m_bps.end_ids[1]:
                        clusters.append(SSandSeqCluster(m_bps.end_ids[1]))
                        clusters[-1].motif_matches(m_bps)

        mes = []
        mes_names = []
        motif_data = []
        motif_keys = ['data', 'name', 'end_name', 'end_id', 'id']
        count = 0
        for c in clusters:
            #if c.end_id != "AG_LL_UU_RR" and  c.end_id != "CU_LL_AG_RR":
            #    continue
            aligned_motifs = []
            for i, m_and_e in enumerate(c.motif_and_ends):
                m, ei = m_and_e.motif, m_and_e.end_index
                m_a = motif_factory.factory.can_align_motif_to_end(m, ei)
                if m_a is None:
                    continue
                m_a = motif_factory.factory.align_motif_to_common_frame(m_a, ei)
                aligned_motifs.append(m_a)
            m_clusters = cluster.cluster_motifs(aligned_motifs)
            clustered_motifs = []
            energies = []
            for j, c_motifs in enumerate(m_clusters):
                clustered_motifs.append(c_motifs.motifs[0])
                pop = float(len(c_motifs.motifs)) / float(len(aligned_motifs))
                energy = -kBT*math.log(pop)
                energies.append(energy)

            me = motif_ensemble.MotifEnsemble()
            me.setup(c.end_id, clustered_motifs, energies)
            mes.append(me)
            mes_names.append(me.id)

            motif = me.members[0].motif
            spl = me.id.split("_")
            motif.name = spl[0][0]+spl[2][1]+"="+spl[0][1]+spl[2][0]

            motif_data.append([motif.to_str(), motif.name, motif.ends[0].name(),
                               me.id, count])
            count += 1



        path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/bp_steps.db"
        sqlite_library.build_sqlite_library(path, mes, mes_names)
        path = settings.RESOURCES_PATH +"/motif_libraries_new/bp_steps.db"
        #sqlite_library.build_sqlite_library(path, motifs, mes_names)
        sqlite_library.build_sqlite_library_2(path, motif_data, motif_keys, 'id')

    def build_motif_state_libraries(self):
        for libname in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            mlib = sqlite_library.MotifSqliteLibrary(libname)
            mlib.load_all()
            motif_states = []
            names = []
            for m in mlib.all():
                ms = m.get_state()
                motif_states.append(ms)
                names.append(ms.name)

            path = settings.RESOURCES_PATH + "/motif_state_libraries/" + libname + ".db"
            sqlite_library.build_sqlite_library(path, motif_states, names)

        for libname in sqlite_library.MotifSSIDSqliteLibrary.get_libnames().keys():
            mlib = sqlite_library.MotifSSIDSqliteLibrary(libname)
            mlib.load_all()
            motif_states = []
            names = []
            for m in mlib.all():
                ms = m.get_state()
                motif_states.append(ms)
                names.append(ms.name)

            path = settings.RESOURCES_PATH + "/motif_state_libraries/ss_" + libname + ".db"
            sqlite_library.build_sqlite_library(path, motif_states, names)

    def build_unique_twoway_library(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        mlib.load_all()
        clusters = cluster.cluster_motifs(mlib.all(), 9.0)
        motif_arrays = []
        motif_array_names = []

        data = []
        keys = ['data', 'name', 'end_name', 'end_id', 'id']

        for i, c in enumerate(clusters):
            lowest = c.motifs[0]
            for m in c.motifs:
                if lowest.score > m.score:
                    lowest = m

            motif_arrays.append(motif.MotifArray(c.motifs))
            motif_array_names.append(lowest.name)

            data.append([lowest.to_str(), lowest.name,
                        lowest.ends[0].name(), lowest.end_ids[0], i])


        path = settings.RESOURCES_PATH +"/motif_libraries_new/unique_twoway.db"
        sqlite_library.build_sqlite_library_2(path, data, keys, 'id')

        #path = settings.RESOURCES_PATH +"/motif_libraries_new/twoway_clusters.db"
        #sqlite_library.build_sqlite_library(path, motif_arrays, motif_array_names)

    def build_ss_and_seq_libraries(self):
        #libnames = ["twoway", "tcontact", "hairpin", "nway"]
        libnames = ["twoway"]

        for libname in libnames:
            mlib = sqlite_library.MotifSqliteLibrary(libname)
            mlib.load_all()

            clusters = []

            mes = []
            mes_names = []

            motifs = []

            kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
            kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

            for m in mlib.all():
                matched = 0
                for c in clusters:
                    if c.motif_matches_end(m, 0):
                        matched = 1

                if not matched:
                    clusters.append(SSandSeqCluster(m.end_ids[0]))
                    clusters[-1].motif_matches_end(m, 0)

            for c in clusters:

                all_motifs = []
                for m_e in c.motif_and_ends:
                    all_motifs.append(m_e.motif)
                if libname != "twoway":
                    energies = [1 for x in all_motifs]
                    me = motif_ensemble.MotifEnsemble()
                    me.setup(c.end_id, all_motifs, energies)
                    mes.append(me)
                    motifs.append(me.members[0].motif)
                    motifs[-1].name = me.id
                    mes_names.append(me.id)
                    continue
                m_clusters = cluster.cluster_motifs(all_motifs)
                clustered_motifs = []
                energies = []
                for j, c_motifs in enumerate(m_clusters):
                    clustered_motifs.append(c_motifs.motifs[0])
                    pop = float(len(c_motifs.motifs)) / float(len(all_motifs))
                    energy = -kBT*math.log(pop)
                    energies.append(energy)

                me = motif_ensemble.MotifEnsemble()
                me.setup(c.end_id, clustered_motifs, energies)
                mes.append(me)
                motifs.append(me.members[0].motif)
                motifs[-1].name = me.id
                mes_names.append(me.id)

            print libname, len(mlib.all()), len(clusters)

            path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/"+libname+".db"
            sqlite_library.build_sqlite_library(path, mes, mes_names)

            path = settings.RESOURCES_PATH +"/motif_libraries_new/ss_"+libname+".db"
            sqlite_library.build_sqlite_library(path, motifs, mes_names)

    def build_motif_ensemble_state_libraries(self):

        for libname in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():

            me_lib = sqlite_library.MotifEnsembleSqliteLibrary(libname)
            me_lib.load_all()

            mses = []
            names = []
            for me in me_lib.all():
                mse = me.get_state()
                mses.append(mse)
                names.append(mse.id)

            path = settings.RESOURCES_PATH +"/motif_state_ensemble_libraries/"+libname+".db"
            sqlite_library.build_sqlite_library(path, mses, names)

#setup_start_motif()

builder = BuildSqliteLibraries()
#builder.build_ideal_helices()
#builder.build_basic_libraries()
#builder.build_helix_ensembles()
#builder.build_ss_and_seq_libraries()
builder.build_unique_twoway_library()
#builder.build_motif_state_libraries()
#builder.build_motif_ensemble_state_libraries()


#mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
#m = mlib.get("HELIX.IDEAL")

#me_lib = sqlite_library.MotifStateEnsembleSqliteLibrary("bp_steps")
#me = me_lib.get("GG_LL_CC_RR")
#print len(me.members)

#print ref_motif.name



#builder.build_ideal_helices()
