import rnamake
import rnamake.settings as settings
import rnamake.motif_library_sqlite as motif_library_sqlite

path = settings.RESOURCES_PATH + "motif_libraries/bp_steps.db"
mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)
mlib.load_all()

motifs = []
names = []
for m in mlib.motifs():
    seq = m.sequence()
    ss = m.secondary_structure()

    seqs = seq.split("&")
    sss = ss.split("&")
    print seq

    name = ""
    for i, s in enumerate(seqs):
        name += s + "_"
        for e in sss[i]:
            if e == "(":
                name += "L"
            elif e == ")":
                name += "R"
            else:
                name += "U"

        if i != len(seqs)-1:
            name += "_"

    if name in names:
        continue
    if name == 'GC_LL_GC_RR':
        print name
    m.name = name
    motifs.append(m)
    names.append(name)

path = settings.RESOURCES_PATH + "motif_libraries"
motif_library_sqlite.build_sqlite_library(path+"/seq_and_ss.db", motifs, names)

path = settings.RESOURCES_PATH + "motif_libraries/seq_and_ss.db"
mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)

m = mlib.get_motif("GA_LL_UC_RR")