import util
import ss_tree

class SS_Chain(object):
    def __init__(self, seq, ss, bounds=[-1,-1]):
        self.seq, self.ss, self.bounds = seq, ss, bounds

    def __repr__(self):
        return "(" + self.seq + ", " +  self.ss + ", " + str(self.bounds) + ")"

    def copy(self):
        return SS_Chain(self.seq, self.ss, self.bounds)

    def flip_ss(self, pos):
        new_ss = ""
        if pos == 0:
            c, f = "(", ")"
        else:
            c, f = ")", "("

        for e in self.ss:
            if e == f:
                new_ss += c
            elif e == c:
                new_ss += f

        self.ss = new_ss

    def to_str(self):
        s = self.seq + ":" + self.ss + ":" + str(self.bounds[0]) + ":" + str(self.bounds[1])
        return s

def str_to_ss_chain(s):
    spl = s.split(":")
    bounds = [ int(spl[2]), int(spl[3])]
    return SS_Chain(spl[0], spl[1], bounds)

class Residue(object):
    def __init__(self, name, num):
        self.name, self.num = name, num

class Chain(object):
    def __init__(self, residues):
        self.residues = residues

    def __repr__(self):
        seq = ""
        for r in self.residues:
            seq += r.name

        return "<Chain: " + seq

    def first(self):
        return self.residues[0]

    def last(self):
        return self.residues[-1]

class Basepair(object):
    def __init__(self, res1, res2):
        self.res1, self.res2 = res1, res2

    def partner(self, r):
        if   r == self.res1:
            return self.res2
        elif r == self.res2:
            return self.res1
        else:
            raise ValueError("call partner with a residue not in basepair")

class Structure(object):
    def __init__(self, chains, basepairs):
        self.chains = chains
        self.basepairs = basepairs
        self.ends = []

    def bp_for_res(self, res):
        for bp in self.basepairs:
            if bp.res1 == res or bp.res2 == res:
                return bp
        return None

    def sequence(self):
        seqs = [x.seq for x in self.ss_chains]
        return "&".join(seqs)

    def secondary_structure(self):
        sss = [x.ss for x in self.ss_chains]
        return "&".join(sss)

    def reorient_ss_and_seq(self, pos1, pos2):
        end = None
        for e in self.ends:
            if e.res1.num == pos1 and e.res2.num == pos2:
                end = e
                break
            if e.res1.num == pos2 and e.res2.num == pos1:
                end = e
                break

        if end is None:
            raise ValueError("there is no ends with that include: " \
                             + str(pos1) + "," + str(pos2))

        all_chains = self.chains[::]
        open_chains = []
        for c in all_chains:
            if c.first() == end.res1 or c.first() == end.res2:
                open_chains.append(c)
                break

        all_chains.remove(open_chains[0])

        if len(open_chains) == 0:
            raise ValueError("could not find chain to start with")

        seen_res = {}
        seen_bp = {}
        saved_bp = None
        structure = ""
        seq = ""
        bounds = [0, 0]
        ss_chains = []
        count = 0
        while len(open_chains) > 0:
            c = open_chains.pop(0)
            bounds = [count+1, count+1]
            for r in c.residues:
                count += 1
                ss = "."
                bp = self.bp_for_res(r)
                saved_bp = None
                if bp is not None:
                    saved_bp = bp
                    partner_res = bp.partner(r)
                    if   bp not in seen_bp and r not in seen_res and \
                         partner_res not in seen_res:
                        seen_res[r] = 1
                        ss = "("
                    elif partner_res in seen_res:
                        if seen_res[partner_res] > 1:
                            ss = "."
                        else:
                            ss = ")"
                            seen_res[r] = 1
                            seen_res[partner_res] += 1

                structure += ss
                seq += r.name

                if saved_bp is not None:
                    seen_bp[saved_bp] = 1

            bounds[1] = count
            ss_chains.append(SS_Chain(seq, structure, bounds))
            structure = ""
            seq = ""


            best_score = 0

            for c in all_chains:
                score = 0
                for r in c.residues:
                    bp = self.bp_for_res(r)
                    if bp in seen_bp:
                        score += 1
                if score > best_score:
                    best_score = score

            best_chains = []
            for c in all_chains:
                score = 0
                for r in c.residues:
                    bp = self.bp_for_res(r)
                    if bp in seen_bp:
                        score += 1
                if score == best_score:
                    best_chains.append(c)

            best_chain = None
            best_score = 10000
            for c in best_chains:
                pos = 1000
                for i, r in enumerate(c.residues):
                    bp = self.bp_for_res(r)
                    if bp in seen_bp:
                        pos = i
                        break
                if pos < best_score:
                    best_score = pos
                    best_chain = c

            if best_chain is None:
                break
            all_chains.remove(best_chain)
            open_chains.append(best_chain)

        self.ss_chains = ss_chains
        return ss_chains

class StructureFactory(object):
    def __init__(self):
        pass



def ss_id(ss_data):
    id = ""
    for i, ss_d in enumerate(ss_data):
        id += ss_d.seq + "_"
        for e in ss_d.ss:
            if   e == "(":
                id += "L"
            elif e == ")":
                id += "R"
            elif e == ".":
                id += "U"
            else:
                raise ValueError("unexpected symbol in dot bracket notation: " + e)
        if i != len(ss_data)-1:
            id += "_"
    return id

def ss_id_to_ss_tree(ss_id):
    ss = ""
    seq = ""
    spl = ss_id.split("_")

    for i in range(0, len(spl)-1, 2):
        seq += spl[i]
        for e in spl[i+1]:
            if   e == "L":
                ss += "("
            elif e == "R":
                ss += ")"
            elif e == "U":
                ss += "."
            else:
                raise ValueError("unexpected symbol in ss_id")

        if i != len(spl)-2:
            seq += "+"
            ss += "+"

    return ss_tree.SS_Tree(ss, seq)


def assign_secondary_structure(motif):
    structure = ""
    seq = ""
    bounds = [0, 0]
    seen_res = {}
    seen_bp = {}
    saved_bp = None
    count = -1
    ss_chains = []

    all_chains = motif.chains()
    open_chains = [all_chains.pop(0)]

    while len(open_chains) > 0:
        c = open_chains.pop(0)
        bounds = [count+1, count+1]
        for r in c.residues:
            count += 1
            ss = ""
            bps = motif.get_basepair(res1=r)
            is_bp = 0
            for bp in bps:
                partner_res = bp.partner(r)
                is_bp = 1
                passes = 0
                saved_bp = None
                if util.wc_bp(bp) and bp.bp_type == "cW-W":
                    passes = 1
                if util.gu_bp(bp) and bp.bp_type == "cW-W":
                    passes = 1

                if passes:
                    saved_bp = bp
                    if   bp not in seen_bp and r not in seen_res and \
                         partner_res not in seen_res:
                        seen_res[r] = 1
                        ss = "("
                    elif partner_res in seen_res:
                        if seen_res[partner_res] > 1:
                            ss = "."
                        else:
                            ss = ")"
                            seen_res[r] = 1
                            seen_res[partner_res] += 1
                            break
                elif r not in seen_res:
                    ss = "."

            if not is_bp:
                ss = "."

            if saved_bp is not None:
                seen_bp[saved_bp] = 1

            structure += ss
            seq += r.name

        bounds[1] = count
        ss_chains.append(SS_Chain(seq, structure, bounds))
        structure = ""
        seq = ""

        best_chain = None
        best_score = 0

        for c in all_chains:
            score = 0
            for r in c.residues:
                bps = motif.get_basepair(res1=r)
                for bp in bps:
                    if bp in seen_bp:
                        score += 1
            if score > best_score:
                best_chain = c
                best_score = score

        if best_chain is None:
            break
        all_chains.remove(best_chain)
        open_chains.append(best_chain)

    return ss_chains