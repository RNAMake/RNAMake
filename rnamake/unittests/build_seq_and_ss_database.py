import rnamake
import rnamake.settings as settings
import rnamake.motif_library_sqlite as motif_library_sqlite
import rnamake.motif_tree as motif_tree
import rnamake.util as util
import rnamake.motif_type as motif_type

def convert_library(mlib):
    mt = motif_tree.MotifTree()
    motifs = []
    names = []
    start = motif_tree.MotifTree().nodes[0].motif
    start_res = start.residues()[1]
    for m in mlib.motifs():
        mt.add_motif(m, end_index=m.end_to_add)
        m = mt.nodes[1].motif
        m._assign_secondary_structure()
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
        if seqs[0] == 'AUAG' or seqs[1] == 'AUAG':
            print "made it"
        m.name = name
        m.mtype = mlib.mtype
        motifs.append(m)
        names.append(name)

    return motifs, names

path = settings.RESOURCES_PATH + "motif_libraries/bp_steps.db"
mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)
mlib.mtype = motif_type.HELIX
mlib.load_all()

motifs1, names1 = convert_library(mlib)

path = settings.RESOURCES_PATH + "motif_libraries/twoway_aligned.db"
mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)
mlib.mtype = motif_type.TWOWAY
mlib.load_all()

motifs2, names2 = convert_library(mlib)

motifs = motifs1 + motifs2
names = names1 + names2


path = settings.RESOURCES_PATH + "motif_libraries"
motif_library_sqlite.build_sqlite_library(path+"/seq_and_ss.db", motifs, names)

path = settings.RESOURCES_PATH + "motif_libraries/seq_and_ss.db"
mlib = motif_library_sqlite.MotifLibrarySqlite(libpath=path)
mlib.mtype = motif_type.HELIX
