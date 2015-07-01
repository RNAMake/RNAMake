import ss_tree
import secondary_structure



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
            ss_chains.append(secondary_structure.SecondaryStructureChain(structure,
                                                                         seq))
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

        return secondary_structure.SecondaryStructure(ss_chains)


class StructureSecondaryFactory(object):
    def __init__(self):
        pass

    def _get_chains(self, seq):
        chains = []
        residues = []

        for i in range(len(seq)):
            if seq[i] != "&":
                r = Residue(seq[i], i)
                residues.append(r)
            else:
                chains.append(Chain(residues))
                residues = []

        if len(residues) > 0:
            chains.append(Chain(residues))

        return chains

    def _get_basepairs(self, sstree, all_res):
        basepairs, ends = [], []
        for n in sstree:
            if n.data.type == ss_tree.SS_Type.SS_BP:
                res1_i = n.data.ss_data[0].bounds[0]
                res2_i = n.data.ss_data[1].bounds[1]
                bp_res = []
                for r in all_res:
                    if r.num == res1_i or r.num == res2_i:
                        bp_res.append(r)
                if len(bp_res) != 2:
                    raise ValueError("incorrect number of resiudes in ss_bp")

                bp = Basepair(bp_res[0], bp_res[1])
                basepairs.append(bp)
                if n.parent is None:
                    ends.append(bp)
                    continue
                children = []
                for c in n.children:
                    if c is None:
                        continue
                    children.append(c)
                if len(children) == 1:
                    if children[0].data.type == ss_tree.SS_Type.SS_SEQ_BREAK:
                        ends.append(bp)
        return basepairs, ends

    def get_structure(self, seq, ss, pos1, pos2):
        sstree = ss_tree.SS_Tree(seq, ss)

        chains = self._get_chains(seq)
        all_res = []
        for c in chains:
            all_res.extend(c.residues)
        basepairs, ends = self._get_basepairs(sstree, all_res)
        struct = Structure(chains, basepairs)
        struct.ends = ends
        ss = struct.reorient_ss_and_seq(pos1, pos2)
        return ss

factory = StructureSecondaryFactory()


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

