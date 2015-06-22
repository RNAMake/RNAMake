import rnamake.motif_library as motif_library
import rnamake.motif_factory as motif_factory
import rnamake.motif_ensemble as motif_ensemble
import rnamake.motif as motif
import rnamake.motif_type as motif_type
import rnamake.util as util
import rnamake.cluster as cluster
import rnamake.settings as settings
import rnamake.sqlite_library as sqlite_library
import math

def test_align(self):
    mlib = sqlite_library.MotifLibrarySqlite(libname="ideal_helices")
    m1 = mlib.get_motif("HELIX.IDEAL.12")
    m2 = m1.copy()

    aligned = motif.get_aligned_motif(m1.ends[1], m2.ends[0], m2)

    m1.to_pdb("test_1.pdb")
    aligned.to_pdb("test_2.pdb")


class SSandSeqCluster(object):

    def __init__(self, end_id):
        self.end_id = end_id
        self.motif_and_ends = []

    def motif_matches(self, m):
        class MotifandEnd(object):
            def __init__(self, m, ei):
                self.motif, self.end_index = m, ei

        for i, end_id in enumerate(m.end_ids):
            if self.end_id == end_id:
                self.motif_and_ends.append(MotifandEnd(m, i))
                return 1
        return 0


class BuildSqliteLibraries(object):
    def __init__(self):
        pass

    def build_ideal_helices(self):
        mlib = motif_library.ideal_helix_lib()
        aligned_motifs = []
        names = []
        path = settings.RESOURCES_PATH +"/motif_libraries_new/ideal_helices.db"
        for m in mlib.motifs():
            aligned_motif = motif_factory.factory.align_motif_to_common_frame(m, 1)
            aligned_motifs.append(m)
            names.append(m.name)
        sqlite_library.build_sqlite_library(path, aligned_motifs, names)

        mlib = sqlite_library.MotifSqliteLibrary(libname="ideal_helices")
        m = mlib.get("HELIX.IDEAL.6")
        s = m.to_str()
        path = settings.RESOURCES_PATH + "/motifs/base.motif"
        f = open(path, "w")
        f.write(s)
        f.close()

    def build_twoways(self):
        mlib = motif_library.unique_twoway_lib()
        h_mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")

        succeses = []
        for m in mlib.motifs():
            if len(m.ends) == 1:
                continue

            for ei in range(len(m.ends)):
                m_added = motif_factory.factory.can_align_motif_to_end(m, ei)
                if m_added is None:
                    continue
                else:
                    succeses.append([m_added, ei])

        aligned_motifs = []
        names = []

        for s in succeses:
            m, ei = s
            m_added = motif_factory.factory.align_motif_to_common_frame(m, ei)
            m_added.name = m_added.name + "-" + m_added.ends[0].name()
            aligned_motifs.append(m_added)
            names.append(m_added.name)


        path = settings.RESOURCES_PATH +"/motif_libraries_new/twoways.db"
        sqlite_library.build_sqlite_library(path, aligned_motifs, names)

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
        motifs = []
        for c in clusters:
            aligned_motifs = []
            for i, m_and_e in enumerate(c.motif_and_ends):
                m, ei = m_and_e.motif, m_and_e.end_index
                m_a = motif_factory.factory.can_align_motif_to_end(m, ei)
                m_a = motif_factory.factory.align_motif_to_common_frame(m_a, ei)
                aligned_motifs.append(m_a)
            clusters = cluster.cluster_motifs(aligned_motifs)
            clustered_motifs = []
            energies = []
            for j, c_motifs in enumerate(clusters):
                clustered_motifs.append(c_motifs.motifs[0])
                pop = float(len(c_motifs.motifs)) / float(len(aligned_motifs))
                energy = -kBT*math.log(pop)
                energies.append(energy)

            me = motif_ensemble.MotifEnsemble()
            me.setup(c.end_id, clustered_motifs, energies)
            mes.append(me)
            motifs.append(me.members[0].motif)
            motifs[-1].name = me.id
            mes_names.append(me.id)

        path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/bp_steps.db"
        sqlite_library.build_sqlite_library(path, mes, mes_names)
        path = settings.RESOURCES_PATH +"/motif_libraries_new/bp_steps.db"
        sqlite_library.build_sqlite_library(path, motifs, mes_names)

    def build_motif_state_libraries(self):
        libnames = ["ideal_helices", "twoways"]
        for libname in libnames:
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









builder = BuildSqliteLibraries()
#builder.build_ideal_helices()
#builder.build_twoways()
#builder.build_helix_ensembles()
builder.build_motif_state_libraries()




#me_lib = sqlite_library.MotifEnsembleSqliteLibrary("bp_steps")
#me = me_lib.get_motif_ensemble("GG_LL_CC_RR")
#print me.id

#print ref_motif.name



#builder.build_ideal_helices()
