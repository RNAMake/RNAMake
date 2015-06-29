import random
import rnamake.resource_manager as rm
import rnamake.motif_tree as motif_tree
import rnamake.ss_tree as ss_tree



class BuildMotifTree(object):
    def __init__(self, lib_names = ["ideal_helices", "twoway"], libs = None):
        if libs is None:
            self.libs = [rm.manager.mlibs[x] for x in lib_names ]
        else:
            self.libs = libs

    def build(self, size = 2):
        mt = motif_tree.MotifTree()
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
                print "fail",m.name
                return mt


        return mt

    def build_specific(self, names):
        size = len(names)
        mt = motif_tree.MotifTree()
        pos = 0
        count = 0
        cpos = 0
        while len(mt) < size:
            if pos == len(self.libs):
                pos = 0

            m = self.libs[pos].get(names[cpos])
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

class BuildSSTree(object):
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

        return ss_tree.SS_Tree(ss, seq)

