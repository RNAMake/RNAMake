import rnamake.motif_ensemble as motif_ensemble
import rnamake.motif as motif
import rnamake.motif_type as motif_type
import rnamake.util as util
import rnamake.cluster as cluster
import rnamake.sqlite_library as sqlite_library
import math
import subprocess
import os
import numpy as np
import shutil

import glob
from rnamake import motif_factory, settings
from collections import defaultdict

bad_keys = "TWOWAY.2GDI.4-X20-X45 TWOWAY.1S72.46-02097-02647 TWOWAY.2GDI.6-Y20-Y45".split()


class BuildSqliteLibraries(object):

    def __init__(self):
        self.mf = motif_factory.MotifFactory()

    def build_ideal_helices(self):
        path = settings.RESOURCES_PATH + "/motif_libraries/ideal_helices.db"
        data = []
        rdata = []
        keys = ['data', 'name', 'end_name', 'end_id', 'id']

        # h_dirs = settings.MOTIF_DIRS + "/helices/"
        h_dirs = "/Users/josephyesselman/projects/RNAMake/rnamake/resources/motifs/helices/"
        ideal_dirs = glob.glob(h_dirs+"HELIX.IDEA*")
        count = 0

        for i_d in ideal_dirs:
            motifs = self.mf.motifs_from_file(i_d, mtype=motif_type.HELIX)
            # f_m the foward ideal helix used in design
            # r_m the reverse ideal helix only used in specific cases
            f_m, r_m = motifs

            data.append([f_m.to_str(), f_m.name,
                         f_m.get_end(0).name, f_m.get_end_id(0), count])

            rdata.append([r_m.to_str(), r_m.name,
                         r_m.get_end(0).name, r_m.get_end_id(0), count])

            count += 1

        sqlite_library.build_sqlite_library(path, data, keys, 'id')
        path = settings.RESOURCES_PATH + "/motif_libraries/ideal_helices_reversed.db"
        sqlite_library.build_sqlite_library(path, rdata, keys, 'id')

    def __correct_ends_for_motif(self, m, t):
        if t == motif_type.HAIRPIN and m.num_ends() != 1:
            return 0
        if t == motif_type.TWOWAY and m.num_ends() != 2:
            return 0
        if t == motif_type.NWAY and m.num_ends() < 3:
            return 0
        if t == motif_type.HELIX and m.num_ends() != 2:
            return 0
        if t == motif_type.TCONTACT and m.num_ends() < 2:
            return 0
        return 1

    def build_basic_libraries(self):

        # types = [motif_type.TWOWAY, motif_type.NWAY, motif_type.HAIRPIN,
        #         motif_type.TCONTACT]

        # FIX on main branch
        base_dirs = "/Users/josephyesselman/projects/RNAMake/rnamake/resources/motifs/"
        types = {
            motif_type.HAIRPIN  : base_dirs + "hairpins",
            motif_type.HELIX    : base_dirs + "helices",
            motif_type.TWOWAY   : base_dirs + "two_ways",
            motif_type.NWAY     : base_dirs + "junctions",
            motif_type.TCONTACT : base_dirs + "tertiary_contacts"
        }

        for t, p in types.iteritems():
            print p
            motif_dirs = glob.glob(p+"/*")

            count = 0
            data = []
            keys = ['data', 'name', 'end_name', 'end_id', 'id']

            for m_dir in motif_dirs:
                if not os.path.isdir(m_dir):
                    continue
                motifs = self.mf.motifs_from_file(m_dir, mtype=t)
                for m in motifs:
                    if m.name + "-" + m.get_end(0).name in bad_keys:
                        continue
                    if not self.__correct_ends_for_motif(m ,t):
                        continue

                    data.append([m.to_str(), m.name, m.get_end(0).name,
                                 m.get_end_id(0), count])
                    count += 1

            path = settings.RESOURCES_PATH +"/motif_libraries/"+\
                   motif_type.type_to_str(t).lower()+".db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def __get_bp_steps(self):
        helix_mlib = sqlite_library.MotifSqliteLibrary("helix")
        helix_mlib.load_all(100)

        all_bp_steps = defaultdict(list)

        seen_name = []
        for m in helix_mlib.all():
            pos = 0
            spl = m.name.split(".")
            if spl[1] == "IDEAL" or spl[1] == "LE":
                continue
            if m.name in seen_name:
                continue
            seen_name.append(m.name)
            for i in range(m.num_basepairs() - 1):
                if m.get_basepair(i).bp_type != "cW-W" or \
                                m.get_basepair(i + 1).bp_type != "cW-W":
                    continue
                bps = [m.get_basepair(i), m.get_basepair(i + 1)]
                res = []
                for bp in bps:
                    bp_res = m.get_bp_res(bp)
                    for r in bp_res:
                        if r not in res:
                            res.append(r)

                if len(res) != 4:
                    continue

                name = m.name + ".BP." + str(pos)
                pos += 1
                bp_steps = self.mf.motifs_from_res(res, bps, m, name, motif_type.HELIX)
                if len(bp_steps) != 2:
                    continue
                for bp_step in bp_steps:
                    if bp_step.num_ends() != 2 or bp_step.num_chains() != 2:
                        break

                    fail = 0
                    for bp in bp_step.iter_basepairs():
                        bp_res = bp_step.get_bp_res(bp)
                        vec1 = bp_res[0].get_coords("C2'") - bp_res[0].get_coords("O4'")
                        vec2 = bp_res[1].get_coords("C2'") - bp_res[1].get_coords("O4'")
                        if np.dot(vec1, vec2) > 2:
                            fail = 1
                            break
                    if fail:
                        break

                    aligned_bp_step = self.mf.align_motif_to_common_frame(bp_step, 0)
                    all_bp_steps[bp_step.get_end_id(0)].append(aligned_bp_step)

        return all_bp_steps

    def __get_renamed_motif(self, m, new_name):
        new_ms = self.mf.motifs_from_rstruc(m, motif_type.HELIX, new_name)
        if new_ms[0].get_end(0).name == m.get_end(0).name:
            new_m = new_ms[0]
        else:
            new_m = new_ms[1]
        return new_m

    def __get_renamed_motifs(self, motifs, bp_count):
        renamed_motifs = []
        for i, m in enumerate(motifs):
            new_m = self.__get_renamed_motif(m, "BP." + str(bp_count) + "." + str(i))
            renamed_motifs.append(new_m)
        return renamed_motifs

    def build_helix_ensembles(self):

        all_bp_steps = self.__get_bp_steps()

        mes_keys = ['data', 'name', 'id']
        mes_data = []

        motif_data = []
        motif_keys = ['data', 'name', 'end_name', 'end_id', 'id']

        unique_data = []

        bp_count = 0
        count = 1

        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

        seen_end_id = []
        for end_id, bp_steps in all_bp_steps.iteritems():
            if end_id in seen_end_id:
                continue
            seen_end_id.append(end_id)
            bp_count += 1
            # best 0.65
            m_clusters = cluster.cluster_motifs(bp_steps, 0.65)

            energies = []
            clustered_motifs = []

            for j, c_motifs in enumerate(m_clusters):
                m = c_motifs.motifs[0]
                clustered_motifs.append(m)

                pop = float(len(c_motifs.motifs)) / float(len(bp_steps))
                energies.append( -kBT*math.log(pop) )

            me = motif_ensemble.MotifEnsemble(self.__get_renamed_motifs(clustered_motifs, bp_count),
                                              energies)
            best_m = self.__get_renamed_motif(me.get_member(0).motif,  "BP."+str(bp_count))

            mes_data.append([me.to_str(), me.end_id, count])

            unique_data.append([best_m.to_str(), best_m.name, best_m.get_end(0).name,
                                best_m.get_end_id(0), count])

            count += 1

            if best_m.get_end_id(1) in seen_end_id:
                continue
            seen_end_id.append(best_m.get_end_id(1))

            other_bp_steps = all_bp_steps[m.get_end_id(1)]
            other_motifs = []
            for m in clustered_motifs:
                for other_bp in other_bp_steps:
                    if m.name == other_bp.name and m.get_end(0).name == other_bp.get_end(1).name:
                        other_motifs.append(other_bp)
                        break

            other_me = motif_ensemble.MotifEnsemble(self.__get_renamed_motifs(other_motifs, bp_steps),
                                                    energies)
            other_m = self.__get_renamed_motif(other_me.get_member(0).motif, "BP."+str(bp_count))

            mes_data.append([other_me.to_str(), other_me.end_id, count])

            unique_data.append([other_m.to_str(), other_m.name, other_m.get_end(0).name,
                                other_m.get_end_id(0), count])
            count += 1

        path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/bp_steps.db"
        sqlite_library.build_sqlite_library(path, mes_data, mes_keys, 'id')
        path = settings.RESOURCES_PATH +"/motif_libraries/bp_steps.db"
        sqlite_library.build_sqlite_library(path, unique_data, motif_keys, 'id')

    def build_new_bp_steps(self):
        mlib = sqlite_library.MotifSqliteLibrary("bp_steps")
        mlib.load_all()

        motifs = []

        motif_data = []
        motif_keys = ['data', 'name', 'end_name', 'end_id', 'id']
        count = 0

        for m in mlib.all():
            spl = m.name.split(".")
            if len(spl) == 1:
                motifs.append(m)

        unique = []
        unique_m = []

        i = 0
        for m in motifs:
            name_spl = m.name.split("=")

            if m.end_ids[0] in unique:
                continue
            if m.end_ids[1] in unique:
                continue

            old_name = m.name
            unique.append( m.end_ids[0])
            if not m.end_ids[1] in unique:
                unique.append( m.end_ids[1])
            #print m.end_ids

            m.name = "BP."+str(i)

            m_a = motif_factory.factory.can_align_motif_to_end(m, 1)
            m_a = motif_factory.factory.align_motif_to_common_frame(m_a, 1)

            motif_data.append([m.to_str(), m.name, m.ends[0].name(),
                               m.end_ids[0], count])

            count += 1

            motif_data.append([m_a.to_str(), m_a.name, m_a.ends[0].name(),
                               m_a.end_ids[0], count])

            unique_m.append(m)
            #print old_name, m.end_ids[0], m_a.end_ids[0]

            count += 1
            i += 1

        print len(unique)
        for m in unique_m:
            print m.end_ids

        path = settings.RESOURCES_PATH +"/motif_libraries_new/new_bp_steps.db"
        sqlite_library.build_sqlite_library(path, motif_data, motif_keys, 'id')

    def build_le_helix_lib(self):


        mt = motif_tree.MotifTree()

        for i in range(20):
            if m is None:
                continue



            data = []
            keys = ['data', 'name', 'end_name', 'end_id', 'id']
            for i, s in enumerate(succeses):
                m, ei = s
                m_added = motif_factory.factory.align_motif_to_common_frame(m, ei)
                #print m_added.name, m_added.ends[0].name(), m_added.end_ids[0]

                data.append([m_added.to_str(), m_added.name,
                            m_added.ends[0].name(), m_added.end_ids[0], i])

            path = settings.RESOURCES_PATH +"/motif_libraries_new/le_helices.db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_motif_state_libraries(self):
        for libname in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            print libname
            data = []
            keys = ['data', 'name', 'end_name', 'end_id', 'id']

            mlib = sqlite_library.MotifSqliteLibrary(libname)
            mlib.load_all()
            motif_states = []
            names = []
            for i, m in enumerate(mlib.all()):
                ms = m.get_state()
                data.append([ms.to_str(), ms.name,
                             ms.end_names[0], ms.end_ids[0], i])



            path = settings.RESOURCES_PATH + "/motif_state_libraries/" + libname + ".db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_unique_twoway_library(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        mlib.load_all()
        clusters = cluster.cluster_motifs(mlib.all(), 9.0)
        motif_arrays = []
        motif_array_names = []

        data = []
        keys = ['data', 'name', 'end_name', 'end_id', 'id']

        mes_keys = ['data', 'name', 'id']
        mes_data = []
        f = open("sim_list_new", "w")

        count = 0
        for i, c in enumerate(clusters):
            lowest = c.motifs[0]
            for m in c.motifs:
                if lowest.score > m.score:
                    lowest = m

            if lowest.name == "TWOWAY.1GID.6" or \
               lowest.name == "TWOWAY.1GID.2" or \
               lowest.name == "TWOWAY.2GDI.4":
                continue
            count += 1

            f.write(lowest.name + "," + lowest.get_end(0).name + " | ")

            for m in c.motifs:
                f.write(m.name + "," + m.get_end(0).name + " | ")
            f.write("\n")

            data.append([lowest.to_str(), lowest.name,
                        lowest.get_end(0).name, lowest.get_end_id(0), count])


        f.close()

        path = settings.RESOURCES_PATH +"/motif_libraries/unique_twoway.db"
        sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_ss_and_seq_libraries(self):
        libnames = ["twoway", "tcontact", "hairpin", "nway"]
        #libnames = ["twoway"]

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

            keys = ['data', 'name', 'id']
            data = []

            for i, c in enumerate(clusters):

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
                    data.append([motifs[-1].to_str(), me.id, i])
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
                data.append([me.to_str(), me.id, i])

            print libname, len(mlib.all()), len(clusters), len(data), len(keys)


            path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/"+libname+".db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_motif_ensemble_state_libraries(self):

        for libname in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
            print libname
            me_lib = sqlite_library.MotifEnsembleSqliteLibrary(libname)
            me_lib.load_all()

            keys = ['data', 'name', 'id']
            data = []

            for i, me in enumerate(me_lib.all()):
                mse = me.get_state()
                #print len(me.members)
                data.append([mse.to_str(), mse.id, i])

            path = settings.RESOURCES_PATH +"/motif_state_ensemble_libraries/"+libname+".db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_trimmed_ideal_helix_library(self):
        for libname in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            if libname != "ideal_helices":
                continue
            data = []
            keys = ['data', 'name', 'end_name', 'end_id', 'id']

            mlib = sqlite_library.MotifSqliteLibrary(libname)
            mlib.load_all()
            motif_states = []
            names = []
            j = 0
            for i, m in enumerate(mlib.all()):
                name = m.name
                spl = name.split(".")
                if len(spl) < 3:
                    continue
                num = int(spl[2])
                if num < 2:
                    continue
                ms = m.get_state()
                data.append([ms.to_str(), ms.name,
                             ms.end_names[0], ms.end_ids[0], j])
                j += 1



            path = settings.RESOURCES_PATH + "/motif_state_libraries/" + libname + "_min.db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

#setup_start_motif()
builder = BuildSqliteLibraries()

#builder.build_ideal_helices()
#builder.build_trimmed_ideal_helix_library()
#builder.build_basic_libraries()
builder.build_helix_ensembles()
#builder.build_new_bp_steps()
#builder.build_ss_and_seq_libraries()
#builder.build_unique_twoway_library()
#builder.build_motif_state_libraries()
#builder.build_motif_ensemble_state_libraries()



#mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
#m = mlib.get("HELIX.IDEAL")

#me_lib = sqlite_library.MotifStateEnsembleSqliteLibrary("bp_steps")
#me = me_lib.get("GG_LL_CC_RR")
#print len(me.members)

#print ref_motif.name



#builder.build_ideal_helices()
