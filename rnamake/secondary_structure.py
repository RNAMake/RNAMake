import util
import ss_tree

class Residue(object):
    def __init__(self, name, dot_bracket, num):
        self.name, self.dot_bracket, self.num = name, dot_bracket, num


class Chain(object):
    def __init__(self, residues=None):
        self.residues = residues
        if self.residues is None:
            self.residues = []

    def __repr__(self):
        seq = ""
        for r in self.residues:
            seq += r.name

        return "<Chain: " + seq

    def first(self):
        return self.residues[0]

    def last(self):
        return self.residues[-1]

    def sequence(self):
        seq = ""
        for r in self.residues:
            seq += r.name
        return seq

    def dot_bracket(self):
        db = ""
        for r in self.residues:
            db += r.dot_bracket
        return db

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


class SecondaryStructure(object):
    def __init__(self, sequence=None, dot_bracket=None, chains=None):
        if sequence is not None and dot_bracket is not None:
            self.chains = self._setup_chains(sequence, dot_bracket)
        else:
            self.chains = chains
            if self.chains is None:
                self.chains = []

    def __repr__(self):
        return "<SecondaryStructure( " + self.sequence() + " " + self.dot_bracket() + ")"

    def _setup_chains(self, sequence, dot_bracket):
        chains = []
        residues = []

        if len(dot_bracket) != len(sequence):
            raise ValueError("sequence and dot bracket are not the same length")

        if dot_bracket[0] != '(' and dot_bracket[0] != '.' and dot_bracket != '&':
            raise ValueError("secondary structure is not valid did you flip seq and ss?")

        count = 0
        for i in range(len(sequence)):
            if sequence[i] != "&" and sequence[i] != "+":
                r = Residue(sequence[i], dot_bracket[i], count)
                residues.append(r)
                count += 1
            else:
                chains.append(Chain(residues))
                residues = []

        if len(residues) > 0:
            chains.append(Chain(residues))


        return chains

    def get_residue(self, num):
        for r in self.residues():
            if r.num == num:
                return r
        raise ValueError("cannot find residue with num " + str(num))

    def residues(self):
        res = []
        for c in self.chains:
            res.extend(c.residues)
        return res

    def sequence(self):
        sequences = [x.sequence() for x in self.chains]
        return "&".join(sequences)

    def dot_bracket(self):
        dot_brackets = [x.dot_bracket() for x in self.chains]
        return "&".join(dot_brackets)

    def id(self):
        id = ""
        for i, chain in enumerate(self.chains):
            id += chain.sequence + "_"
            for e in chain.dot_bracket:
                if   e == "(":
                    id += "L"
                elif e == ")":
                    id += "R"
                elif e == ".":
                    id += "U"
                else:
                    raise ValueError("unexpected symbol in dot bracket notation: " + e)
            if i != len(self.chains)-1:
                id += "_"
        return id


def str_to_ss_chain(s):
    spl = s.split(":")
    bounds = [ int(spl[2]), int(spl[3])]
    return SS_Chain(spl[0], spl[1], bounds)

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