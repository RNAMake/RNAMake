import rnamake
import rnamake.settings as settings
import rnamake.motif_library_sqlite as motif_library_sqlite
import rnamake.motif_tree as motif_tree
import rnamake.util as util
import rnamake.motif_type as motif_type

path = settings.RESOURCES_PATH + "motif_libraries/bp_steps.db"
mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)
mlib.load_all()

motifs = []
names = []

start = motif_tree.MotifTree().nodes[0].motif
start_res = start.residues()[1]


mt = motif_tree.MotifTree()


for m in mlib.motifs():
    mt.add_motif(m)
    m = mt.nodes[1].motif
    mt.remove_node(mt.last_node)

    chains = m.chains()
    closest = None
    best = 10000
    for c in chains:
        c1 = util.center(c.first().atoms)
        c2 = util.center(start_res.atoms)
        dist = util.distance(c1, c2)
        if dist < best:
            best = dist
            closest = c


    updated_chains = [closest]
    for c in chains:
        if c != closest:
            updated_chains.append(c)

    m.structure.chains = updated_chains
    m.structure._cache_coords()

    seq = m.sequence()
    ss = m.secondary_structure()


    seqs = seq.split("&")
    sss = ss.split("&")
    print seqs

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
    m.name = name
    m.mtype =  motif_type.HELIX
    motifs.append(m)
    names.append(name)

path = settings.RESOURCES_PATH + "motif_libraries"
motif_library_sqlite.build_sqlite_library(path+"/seq_and_ss.db", motifs, names)

path = settings.RESOURCES_PATH + "motif_libraries/seq_and_ss.db"
mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)
mlib.mtype = motif_type.HELIX
