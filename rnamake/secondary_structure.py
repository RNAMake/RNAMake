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
            else:
                new_ss += e

        self.ss = new_ss

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
    for c in motif.chains():
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
            seq += r.rtype.name[0]
        bounds[1] = count
        ss_chains.append(SS_Chain(seq, structure, bounds))
        structure = ""
        seq = ""
    return ss_chains