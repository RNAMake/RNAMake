import util
import ss_tree

class SecondaryStructureChain(object):
    def __init__(self, dot_bracket="", sequence=""):
        self.dot_bracket, self.sequence = dot_bracket, sequence

class SecondaryStructure(object):
    def __init__(self, chains=[]):
        self.chains = chains

    def __repr__(self):
        return "<SecondaryStructure( " + self.sequence() + " " + self.dot_bracket() + ")"

    def sequence(self):
        sequences = [x.sequence for x in self.chains]
        return "&".join(sequences)

    def dot_bracket(self):
        dot_brackets = [x.dot_bracket for x in self.chains]
        return "&".join(dot_brackets)




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