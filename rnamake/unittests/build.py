import random

from rnamake import motif_graph, motif_tree, motif_type, sqlite_library


def fill_basepairs_in_ss(ss):
    pairs = ["AU", "UA", "GC", "CG"]
    for bp in ss.basepairs:
        if bp.res1.name == "N" and bp.res2.name == "N":
            p = random.choice(pairs)
            bp.res1.name = p[0]
            bp.res2.name = p[1]

    for m in ss.motifs:
        for i, end in enumerate(m.ends):
            m.end_ids[i] = secondary_structure.assign_end_id_new(m, end)


class BuildMotifTree(object):
    def __init__(self, rm, lib_names = ["ideal_helices", "unique_twoway"], libs = None):
        self.rm = rm
        if libs is None:
            self.libs = [sqlite_library.MotifSqliteLibrary(x) for x in lib_names ]
        else:
            self.libs = libs

    def build(self, size=2):
        mt = motif_tree.MotifTree(self.rm)
        pos = 0
        count = 0
        while len(mt) < size:
            if pos == len(self.libs):
                pos = 0

            m = self.libs[pos].get_random()
            i = mt.add_motif(m)
            count += 1
            if count > 100:
                break
            if i != -1:
                pos += 1
            else:
                #print "fail",m.name
                return mt
        return mt

    def build_no_ideal_helices(self, size=2):
        mt = self.build(size)
        ss = mt.designable_secondary_structure()
        fill_basepairs_in_ss(ss)
        conn = ss.motif_topology_from_end(ss.ends[0])
        mtt = motif_tree_topology.MotifTreeTopology(conn)
        return motif_tree.motif_tree_from_topology(mtt)

    def build_specific(self, names):
        size = len(names)
        mt = motif_tree.MotifTree()
        pos = 0
        count = 0
        cpos = 0
        while len(mt) < size:
            if pos == len(self.libs):
                pos = 0

            m = self.libs[pos].get(name=names[cpos])
            i = mt.add_motif(m)
            count += 1
            if count > 100:
                break
            if i != -1:
                pos += 1
                cpos += 1
            else:
                print "fail",m.name
                return mt

        return mt


class BuildMotifGraph(object):
    def __init__(self, lib_names = ["ideal_helices", "unique_twoway"], libs = None):
        if libs is None:
            self.libs = [rm.manager.mlibs[x] for x in lib_names ]
        else:
            self.libs = libs

    def build(self, size=2):
        mg = motif_graph.MotifGraph()
        pos = 0
        count = 0
        while len(mg) < size:
            if pos == len(self.libs):
                pos = 0

            m = self.libs[pos].get_random()
            i = mg.add_motif(m)
            count += 1
            if count > 100:
                break
            if i != -1:
                pos += 1
            else:
                #print "fail",m.name
                return mg
        return mg


class BuildSecondaryStructure(object):
    def __init__(self):
        pass

    def build_helix(self, size=10):
        pairs = "AU,UA,GC,CG".split(",")
        s1, s2, ss1, ss2 = "", "", "", ""
        for i in range(size):
            pair = random.choice(pairs)
            s1 += pair[0]
            s2 = pair[1] + s2
            ss1 += "("
            ss2 += ")"

        seq = s1  + "+" + s2
        ss  = ss1 + "+" + ss2

        return ssfactory.factory.pose(seq, ss)


























