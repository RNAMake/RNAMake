import util
import ss_tree
import Queue

class Residue(object):
    def __init__(self, name, dot_bracket, num):
        self.name, self.dot_bracket, self.num = name, dot_bracket, num

    def to_str(self):
        return self.name + "," + self.dot_bracket + "," + str(self.num)

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

    def to_str(self):
        s = ""
        for r in self.residues:
            s += r.to_str() + ";"
        return s


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


class SecondaryStructureMotif(object):
    def __init__(self, type, ends, chains):
        self.type, self.ends, self.chains = type, ends, chains
        self.basepairs = []
        self.end_ids = []

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

    def get_bp(self, res1, res2=None):
        for bp in self.basepairs:
            if bp.res1 != res1 and bp.res2 != res1:
                continue
            if res2 is not None and bp.res1 != res2 and bp.res2 != res2:
                continue
            return bp
        return None

    def to_str(self):
        s = self.type + ";"
        for c in self.chains:
            for r in c.residues:
                s += str(r.num) + " "
            s += ","
        s += ";"
        for bp in self.basepairs:
            s += str(bp.res1.num) + " " + str(bp.res2.num) + ","
        s += ";"
        for end in self.ends:
            s += str(bp.res1.num) + " " + str(bp.res2.num) + ","
        return s


class SecondaryStructure(SecondaryStructureMotif):
    def __init__(self, sequence=None, dot_bracket=None, chains=None):
        if sequence is not None and dot_bracket is not None:
            self.chains = self._setup_chains(sequence, dot_bracket)
        else:
            self.chains = chains
            if self.chains is None:
                self.chains = []

        self.elements = {}
        self.basepairs = []
        self.ends = []

    def __repr__(self):
        return "<SecondaryStructure( " + self.sequence() + " " + self.dot_bracket() + " )"

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

    def motif_topology_from_end(self, end):
        elements = []
        for eles in self.elements.itervalues():
            elements.extend(eles)

        seen_e = {}
        seen_ends = {}
        current_e = None
        for e in elements:
            for bp in e.ends:
                if bp == end:
                    current_e = e
                    break
            if current_e is not None:
                break

        if current_e is None:
            raise ValueError("could not find element with current end")

        current_end = end
        count = 0
        connectivity = []
        queue = Queue.Queue()
        queue.put([current_e, current_end, -1])
        pos = 0
        indexed_e = {}
        while not queue.empty():
            seen_e[current_e] = pos
            indexed_e[pos] = current_e
            current_e, current_end, parent_pos = queue.get()
            ss_id = assign_end_id(current_e, current_end)

            parent_id = ""
            if pos != 0:
                parent = indexed_e[parent_pos]
                parent_id = assign_end_id(parent, current_end)

            connectivity.append([ss_id, parent_id, parent_pos])

            new_ends = []
            for end in current_e.ends:
                if end != current_end:
                    new_ends.append(end)

            for e in elements:
                if e in seen_e:
                    continue
                if e == current_e:
                    continue
                for end in e.ends:
                    if end in new_ends:
                        queue.put([e, end, pos])

            pos += 1

        return connectivity

    def to_str(self):
        s = ""
        for c in self.chains:
            s += c.to_str() + ":"
        s += "@"
        for bp in self.basepairs:
            s += str(bp.res1.num) + " " + str(bp.res2.num) + ","
        s += "@"
        for end in self.ends:
            s += str(end.res1.num) + " " + str(end.res2.num) + ","
        s += "@"
        for eles in self.elements.itervalues():
            for e in eles:
                s += e.to_str() + "|"
        return s


def str_to_residue(s):
    spl = s.split(",")
    return Residue(spl[0], spl[1], int(spl[2]))


def str_to_chain(s):
    """
    creates a chain from string generated from chain.to_str()
    """
    spl = s.split(";")
    c = Chain()
    for r_str in spl[:-1]:
        r = str_to_residue(r_str)
        c.residues.append(r)
    return c


def str_to_secondary_structure(s):
    spl = s.split("@")
    c_spl = spl[0].split(":")
    chains = []
    for chain_str in c_spl[:-1]:
        c = str_to_chain(chain_str)
        chains.append(c)
    ss = SecondaryStructure(chains=chains)
    bp_spl = spl[1].split(",")
    basepairs = []
    for bp_str in bp_spl[:-1]:
        res_nums = bp_str.split()
        res1 = ss.get_residue(int(res_nums[0]))
        res2 = ss.get_residue(int(res_nums[1]))
        bp = Basepair(res1, res2)
        basepairs.append(bp)
    ss.basepairs = basepairs
    end_spl = spl[2].split(",")
    ends = []
    for end_str in end_spl[:-1]:
        res_nums = end_str.split()
        res1 = ss.get_residue(int(res_nums[0]))
        res2 = ss.get_residue(int(res_nums[1]))
        bp = ss.get_bp(res1, res2)
        ends.append(bp)
    motif_spl = spl[3].split("|")
    for motif_str in motif_spl[:-1]:
        motif_info = motif_str.split(";")
        type = motif_info[0]
        chain_spl = motif_info[1].split(",")
        chains = []
        for res_list in chain_spl[:-1]:
            res = res_list.split()
            residues = []
            for r_num in res:
                r = ss.get_residue(int(r_num))
                residues.append(r)
            chains.append(Chain(residues))
        basepairs = []
        bp_spl = motif_str.split(",")
        for bp_str in bp_spl[:-1]:
            res_num = bp_str.split()
            res1 = ss.get_residue(int(res_nums[0]))
            res2 = ss.get_residue(int(res_nums[1]))
            bp = ss.get_bp(res1, res2)
            basepairs.append(bp)
        motif = SecondaryStructureMotif(type, ends, chains)
        motif.basepairs = basepairs
        end_spl = motif_str.split(",")
        ends = []
        for bp_str in bp_spl[:-1]:
            res_num = bp_str.split()
            res1 = ss.get_residue(int(res_nums[0]))
            res2 = ss.get_residue(int(res_nums[1]))
            bp = ss.get_bp(res1, res2)
            ends.append(bp)
        if type not in ss.elements:
            ss.elements[type] = []
        ss.elements[type].append(motif)

    return ss


def assign_secondary_structure(motif):
    bounds = [0, 0]
    seen_res = {}
    seen_bp = {}
    saved_bp = None
    count = -1
    ss_chains = []

    all_chains = motif.chains()[::]
    open_chains = [all_chains.pop(0)]

    while len(open_chains) > 0:
        c = open_chains.pop(0)
        ss_res = []

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

            ss_res.append(Residue(r.name, ss, r.num))
        ss_chains.append(Chain(ss_res))

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

    return SecondaryStructure(chains=ss_chains)


def assign_end_id(ss, end):
    if end not in ss.ends:
        raise ValueError("supplied an end that is not in current ss element")


    all_chains = ss.chains[::]
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
            dot_bracket = "."
            bp = ss.get_bp(r)
            saved_bp = None
            if bp is not None:
                saved_bp = bp
                partner_res = bp.partner(r)
                if   bp not in seen_bp and r not in seen_res and \
                                partner_res not in seen_res:
                    seen_res[r] = 1
                    dot_bracket = "("
                elif partner_res in seen_res:
                    if seen_res[partner_res] > 1:
                        dot_bracket = "."
                    else:
                        dot_bracket = ")"
                        seen_res[r] = 1
                        seen_res[partner_res] += 1

            structure += dot_bracket
            seq += r.name

            if saved_bp is not None:
                seen_bp[saved_bp] = 1

        bounds[1] = count
        ss_chains.append([seq, structure])
        structure = ""
        seq = ""


        best_score = 0

        for c in all_chains:
            score = 0
            for r in c.residues:
                bp = ss.get_bp(r)
                if bp is not None and bp in seen_bp:
                    score += 1
            if score > best_score:
                best_score = score

        best_chains = []
        for c in all_chains:
            score = 0
            for r in c.residues:
                bp = ss.get_bp(r)
                if bp is not None and bp in seen_bp:
                    score += 1
            if score == best_score:
                best_chains.append(c)

        best_chain = None
        best_score = 10000
        for c in best_chains:
            pos = 1000
            for i, r in enumerate(c.residues):
                bp = ss.get_bp(r)
                if bp is not None and bp in seen_bp:
                    pos = i
                    break
            if pos < best_score:
                best_score = pos
                best_chain = c

        if best_chain is None:
            break
        all_chains.remove(best_chain)
        open_chains.append(best_chain)

    ss_id = ""
    for i, chain in enumerate(ss_chains):
        ss_id += chain[0] + "_"
        for e in chain[1]:
            if   e == "(":
                ss_id += "L"
            elif e == ")":
                ss_id += "R"
            elif e == ".":
                ss_id += "U"
            else:
                raise ValueError("unexpected symbol in dot bracket notation: " + e)
        if i != len(ss_chains)-1:
            ss_id += "_"
    return ss_id