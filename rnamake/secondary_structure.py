import util
import ss_tree
import Queue
import uuid
import tree
import motif_type

class Residue(object):
    def __init__(self, name, dot_bracket, num, chain_id, uuid, i_code=""):
        self.name, self.dot_bracket, self.num = name, dot_bracket, num
        self.uuid = uuid
        self.chain_id, self.i_code = chain_id, i_code

    def to_str(self):
        return self.name + "," + self.dot_bracket + "," + str(self.num) + "," + \
               str(self.chain_id) + "," + str(self.i_code)

    def copy(self):
        return Residue(self.name, self.dot_bracket, self.num, self.chain_id, self.uuid, self.i_code)


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

    def copy(self):
        residues = []
        for r in self.residues:
            residues.append(r.copy())
        return Chain(residues)


class Basepair(object):
    def __init__(self, res1, res2, bp_uuid=None):
        self.res1, self.res2 = res1, res2
        self.uuid = bp_uuid
        if self.uuid is None:
            self.uuid = uuid.uuid1()

    def __repr__(self):
        return "<Basepair("+self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code) +\
            "-" + self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code) + ")>"

    def name(self):
        str1 = self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code)
        str2 = self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

        if str1 < str2:
            return str1+"-"+str2
        else:
            return str2+"-"+str1

    def partner(self, r):
        if   r == self.res1:
            return self.res2
        elif r == self.res2:
            return self.res1
        else:
            raise ValueError("call partner with a residue not in basepair")


class Structure(object):
    def __init__(self, chains=None, sequence="", dot_bracket=""):
        self.chains = []
        if chains is not None:
            self.chains = chains

        if len(sequence) != 0 and len(dot_bracket) != 0:
           self.chains = self._setup_chains(sequence, dot_bracket)

    def _setup_chains(self, sequence, dot_bracket):
        chains = []
        residues = []

        if len(dot_bracket) != len(sequence):
            raise ValueError("sequence and dot bracket are not the same length")

        if dot_bracket[0] != '(' and dot_bracket[0] != '.' and dot_bracket != '&':
            raise ValueError("secondary structure is not valid did you flip seq and ss?")

        count = 1
        chains_ids = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K"]
        chain_i = 0
        for i in range(len(sequence)):
            if sequence[i] != "&" and sequence[i] != "+" and sequence[i] != "-":
                r = Residue(sequence[i], dot_bracket[i], count,
                            chains_ids[chain_i], uuid.uuid1())
                residues.append(r)
                count += 1
            else:
                chain_i += 1
                chains.append(Chain(residues))
                residues = []

        if len(residues) > 0:
            chains.append(Chain(residues))

        return chains

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

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        found = []
        for r in self.residues():
            if num is not None and num != r.num:
                continue
            if i_code is not None and i_code != r.i_code:
                continue
            if chain_id is not None and chain_id != r.chain_id:
                continue
            if uuid is not None and uuid != r.uuid:
                continue
            found.append(r)

        if len(found) == 0:
            return None

        if len(found) > 1:
            #print "num,chain_id,icode=",num,chain_id,i_code,
            raise ValueError(
                "found multiple residues in get_residue(), narrow " +
                "your search")

        return found[0]

    def copy(self):
        new_chains = [ c.copy() for c in self.chains]
        return Structure(chains=new_chains)

    def to_str(self):
        s = ""
        for c in self.chains:
            s += c.to_str() + "|"
        return s


class RNAStructure(object):
    def __init__(self, structure=None, basepairs=None, ends=None, name="assembled",
                 path="assembled", score=0, end_ids=None):
        self.structure = structure
        if self.structure is None:
            self.structure = Structure()
        self.basepairs =basepairs
        if self.basepairs is None:
            self.basepairs = []
        self.name = name
        self.path = path
        self.score = score
        self.ends = ends
        if self.ends is None:
            self.ends = []
        self.end_ids = end_ids
        if self.end_ids is None:
            self.end_ids = []

    def __repr__(self):
        return "<secondary_structure.RNAStructure( " + self.sequence() + " " + self.dot_bracket() + " )"

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        return self.structure.get_residue(num, chain_id, i_code, uuid)

    def get_basepair(self, res1=None, res2=None, uuid=None, name=None):
        for bp in self.basepairs:
            if res1 is not None and (bp.res1 != res1 and bp.res2 != res1):
                continue
            if res2 is not None and (bp.res1 != res2 and bp.res2 != res2):
                continue
            if uuid is not None and bp.uuid != uuid:
                continue
            if name is not None and bp.name() != name:
                continue
            return bp
        return None

    def sequence(self):
        return self.structure.sequence()

    def dot_bracket(self):
        return self.structure.dot_bracket()

    def replace_sequence(self, seq):
        spl = seq.split("&")
        seq2 = "".join(spl)
        for i, r in enumerate(self.structure.residues()):
            r.name = seq2[i]

    def residues(self):
        return self.structure.residues()

    def chains(self):
        return self.structure.chains

    def copy(self):
        n_ss = self.structure.copy()
        basepairs, ends = [], []
        for bp in self.basepairs:
            new_bp = Basepair(n_ss.get_residue(uuid=bp.res1.uuid),
                              n_ss.get_residue(uuid=bp.res2.uuid),
                              bp.uuid)
            basepairs.append(new_bp)
        for end in self.ends:
            i = self.basepairs.index(end)
            ends.append(basepairs[i])

        return RNAStructure(n_ss, basepairs, ends, self.name, self.path,
                            self.score, self.end_ids[::])


class Motif(RNAStructure):
    def __init__(self, structure=None, basepairs=None, ends=None,
                 name="assembled", path="assembled", mtype=motif_type.UNKNOWN,
                 score=0, end_ids=None,  id=uuid.uuid1(), r_struct=None):
        self.structure = structure
        if self.structure is None:
            self.structure = Structure()
        self.basepairs = basepairs
        if self.basepairs is None:
            self.basepairs = []
        self.name = name
        self.path = path
        self.score = score
        self.ends = ends
        if self.ends is None:
            self.ends = []
        self.end_ids = end_ids
        if self.end_ids is None:
            self.end_ids = []
        self.mtype = mtype
        self.id = id
        if r_struct is not None:
            self.__dict__.update(r_struct.__dict__)

    def __repr__(self):
        return "<secondary_structure.Motif( " + self.sequence() + " " + self.dot_bracket() + " )"

    def copy(self):
        n_ss = self.structure.copy()
        basepairs, ends = [], []
        for bp in self.basepairs:
            new_bp = Basepair(n_ss.get_residue(uuid=bp.res1.uuid),
                              n_ss.get_residue(uuid=bp.res2.uuid),
                              bp.uuid)
            basepairs.append(new_bp)
        for end in self.ends:
            i = self.basepairs.index(end)
            ends.append(basepairs[i])

        return Motif(n_ss, basepairs, ends, self.name, self.path, self.mtype,
                     self.score, self.end_ids[::], self.id)

    def copy_w_res(self, res, bps):
        chains = []
        m = Motif()
        for c in self.structure.chains:
            new_res = []
            for r in c.residues:
                new_r = res[r.uuid]
                new_res.append(new_r)
            chains.append(Chain(new_res))

        m.structure.chains = chains
        m.end_ids = self.end_ids[::]
        new_bps = []
        new_ends = []
        for bp in self.basepairs:
            new_bps.append(bps[bp.uuid])
        for end in self.ends:
            new_ends.append(bps[bp.uuid])
        m.basepairs = new_bps
        m.ends = new_ends
        return m

    def to_str(self):
        s = str(self.mtype) + "!" + self.name + "!" + self.path + "!" + self.structure.to_str() + "!"
        res = self.residues()
        for bp in self.basepairs:
            s += str(res.index(bp.res1)) + " " + str(res.index(bp.res2)) + "@"
        s += "!"
        for end in self.ends:
            s += str(self.basepairs.index(end)) + " "
        s += "!"
        for ei in self.end_ids:
            s += ei + " "
        return s


class Pose(RNAStructure):
    def __init__(self, structure=None, basepairs=None, ends=None,
                 name="assembled", path="assembled", score=0, end_ids=None,
                 r_struct=None):
        self.structure = structure
        if self.structure is None:
            self.structure = Structure()
        self.basepairs = basepairs
        if self.basepairs is None:
            self.basepairs = []
        self.name = name
        self.path = path
        self.score = score
        self.ends = ends
        if self.ends is None:
            self.ends = []
        self.end_ids = end_ids
        if self.end_ids is None:
            self.end_ids = []
        self.motifs = []

        if r_struct is not None:
            self.__dict__.update(r_struct.__dict__)

    def __repr__(self):
        return "<secondary_structure.Pose( " + self.sequence() + " " + self.dot_bracket() + " )"

    def motif(self, m_id):
        for m in self.motifs:
            if m.id == m_id:
                return m
        return None

    def replace_sequence(self, seq):
        super(self.__class__, self).replace_sequence(seq)

        for m in self.motifs:
            for i, end in enumerate(m.ends):
                m.end_ids[i] = assign_end_id_new(m, end)

    def copy(self):
        c_rna_struct = super(self.__class__, self).copy()
        c_p = Pose(r_struct=c_rna_struct)
        res = {r.uuid  : r for r in self.residues()}
        bps = {bp.uuid : bp for bp in self.basepairs}
        new_motifs = []
        for m in self.motifs:
            m_copy = m.copy_w_res(res, bps)
            new_motifs.append(m_copy)
        c_p.motifs = new_motifs
        return c_p

    def to_str(self):
        s = self.name + "#" + self.path + "#" + self.structure.to_str() + "#"
        res = self.residues()
        for bp in self.basepairs:
            s += str(res.index(bp.res1)) + " " + str(res.index(bp.res2)) + "@"
        s += "#"
        for end in self.ends:
            s += str(self.basepairs.index(end)) + " "
        s += "#"
        for ei in self.end_ids:
            s += ei + " "
        s += "#"
        for m in self.motifs:
            s += m.to_str() + "#"
        return s

    def build_helices(self):
        pass


def str_to_residue(s):
    spl = s.split(",")
    return Residue(spl[0], spl[1], int(spl[2]), spl[3], uuid.uuid1(), spl[4])


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


def str_to_structure(s):
    spl = s.split("|")
    chains = []
    for c_str in spl[:-1]:
        c = str_to_chain(c_str)
        chains.append(c)
    return Structure(chains)


def str_to_motif(s):
    spl = s.split("!")
    m = Motif()
    m.mtype = int(spl[0])
    m.name = spl[1]
    m.path = spl[2]
    m.structure = str_to_structure(spl[3])
    res = m.residues()
    for bp_str in spl[4].split("@")[:-1]:
        res_is = bp_str.split()
        r1 = res[int(res_is[0])]
        r2 = res[int(res_is[1])]
        m.basepairs.append(Basepair(r1, r2))
    for end_i in spl[5].split():
        m.ends.append(m.basepairs[int(end_i)])
    m.end_ids = spl[6].split()

    return m


def str_to_pose(s):
    spl = s.split("#")
    p = Pose()
    p.name = spl[0]
    p.path = spl[1]
    p.structure = str_to_structure(spl[2])
    res = p.residues()
    for bp_str in spl[3].split("@")[:-1]:
        res_is = bp_str.split()
        r1 = res[int(res_is[0])]
        r2 = res[int(res_is[1])]
        p.basepairs.append(Basepair(r1, r2))
    for end_i in spl[4].split():
        p.ends.append(p.basepairs[int(end_i)])
    p.end_ids = spl[5].split()

    motifs = []
    for str in spl[6:-1]:
        spl2 = str.split("!")
        m = Motif()
        m.mtype = int(spl2[0])
        m.name = spl2[1]
        m.path = spl2[2]
        m.structure = str_to_structure(spl2[3])
        chains = []
        for c in m.chains():
            res = []
            for r in c.residues:
                r_new = p.get_residue(r.num, r.chain_id)
                res.append(r_new)
            chains.append(Chain(res))
        m.structure.chains = chains
        res = m.residues()
        for bp_str in spl2[4].split("@")[:-1]:
            res_is = bp_str.split()
            r1 = res[int(res_is[0])]
            r2 = res[int(res_is[1])]
            for bp in p.basepairs:
                if bp.res1 == r1 and bp.res2 == r2:
                    m.basepairs.append(bp)
        for end_i in spl2[5].split():
            m.ends.append(m.basepairs[int(end_i)])
        m.end_ids = spl2[6].split()
        motifs.append(m)
    p.motifs = motifs
    return p


def assign_end_id_new(ss, end):
    if end not in ss.ends:
        raise ValueError("supplied an end that is not in current ss element")

    all_chains = ss.structure.chains[::]
    open_chains = []
    for c in all_chains:
        if c.first() == end.res1 or c.first() == end.res2:
            open_chains.append(c)
            break

    if len(open_chains) == 0:
        raise ValueError("could not find chain to start with")

    all_chains.remove(open_chains[0])

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
            bp = ss.get_basepair(r)
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


        best_score = -1

        for c in all_chains:
            score = 0
            for r in c.residues:
                bp = ss.get_basepair(r)
                if bp is not None and bp in seen_bp:
                    score += 1
            if score > best_score:
                best_score = score

        best_chains = []
        for c in all_chains:
            score = 0
            for r in c.residues:
                bp = ss.get_basepair(r)
                if bp is not None and bp in seen_bp:
                    score += 1
            if score == best_score:
                best_chains.append(c)

        best_chain = None
        best_score = 10000
        for c in best_chains:
            pos = 1000
            for i, r in enumerate(c.residues):
                bp = ss.get_basepair(r)
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

