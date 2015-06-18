import rnamake.motif_library as motif_library
import rnamake.motif_library_sqlite as motif_library_sqlite
import rnamake.motif_factory as motif_factory
import rnamake.motif as motif
import rnamake.util as util
import rnamake.settings as settings

def test_align(self):
    mlib = motif_library_sqlite.MotifLibrarySqlite(libname="ideal_helices")
    m1 = mlib.get_motif("HELIX.IDEAL.12")
    m2 = m1.copy()

    aligned = motif.get_aligned_motif(m1.ends[1], m2.ends[0], m2)

    m1.to_pdb("test_1.pdb")
    aligned.to_pdb("test_2.pdb")


class BuildSqliteLibraries(object):
    def __init__(self):
        pass

    def build_ideal_helices(self):
        mlib = motif_library.ideal_helix_lib()
        aligned_motifs = []
        path = settings.RESOURCES_PATH +"/motif_libraries_new/ideal_helices.db"
        for m in mlib.motifs():
            aligned_motif = motif_factory.factory.align_motif_to_common_frame(m, 1)
            aligned_motifs.append(m)
        motif_library_sqlite.build_sqlite_library(path, aligned_motifs)

    def build_twoways(self):
        mlib = motif_library.unique_twoway_lib()
        h_mlib = motif_library_sqlite.MotifLibrarySqlite(libname="ideal_helices")

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


        path = settings.RESOURCES_PATH +"/motif_libraries_new/twoways.db"
        motif_library_sqlite.build_sqlite_library(path, aligned_motifs)


builder = BuildSqliteLibraries()
builder.build_twoways()


#print ref_motif.name



#builder.build_ideal_helices()
