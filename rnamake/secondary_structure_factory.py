import ss_tree
import secondary_structure
import util


class MotiftoSecondaryStructure(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.chains = []
        self.open_chains = []
        self.seen_res = {}
        self.seen_bp = {}

    def _get_next_chain(self, motif):
        best_score = -1

        for c in self.chains:
            score = 0
            for r in c.residues:
                bps = motif.get_basepair(res1=r)
                for bp in bps:
                    if bp in self.seen_bp:
                        score += 1
            if score > best_score:
                best_score = score

        best_chains = []
        for c in self.chains:
            score = 0
            for r in c.residues:
                bps = motif.get_basepair(res1=r)
                for bp in bps:
                    if bp in self.seen_bp:
                        score += 1
            if score == best_score:
                best_chains.append(c)

        best_chain = None
        best_score = 10000
        for c in best_chains:
            pos = 1000
            for i, r in enumerate(c.residues):
                bps = motif.get_basepair(res1=r)
                for bp in bps:
                    if bp in self.seen_bp:
                        pos = i
                        break
            if pos < best_score:
                best_score = pos
                best_chain = c

        return best_chain

    def _setup_basepairs_and_ends(self, ss, motif):
        ss_bps = []
        for bp in self.seen_bp.keys():
            res1 = ss.get_residue(uuid=bp.res1.uuid)
            res2 = ss.get_residue(uuid=bp.res2.uuid)
            ss_bps.append(secondary_structure.Basepair(res1, res2))
        ss.basepairs = ss_bps
        ss_ends = []

        for end in motif.ends:
            res1 = ss.get_residue(uuid=end.res1.uuid)
            res2 = ss.get_residue(uuid=end.res2.uuid)
            bp = ss.get_bp(res1, res2)
            if bp is None:
                raise ValueError("did not properly find end in generating ss")
            ss_ends.append(bp)
        ss.ends = ss_ends

    def to_secondary_structure(self, motif):
        saved_bp = None
        ss_chains = []

        self.chains = motif.chains()[::]
        self.open_chains = [self.chains.pop(0)]

        while len(self.open_chains) > 0:
            c = self.open_chains.pop(0)
            ss_res = []

            for r in c.residues:
                ss = "."
                bps = motif.get_basepair(res1=r)
                is_bp = 0
                for bp in bps:
                    partner_res = bp.partner(r)
                    passes = 0
                    saved_bp = None
                    if util.wc_bp(bp) and bp.bp_type == "cW-W":
                        passes = 1
                    if util.gu_bp(bp) and bp.bp_type == "cW-W":
                        passes = 1

                    if passes:
                        saved_bp = bp
                        if   bp not in self.seen_bp and \
                              r not in self.seen_res and \
                              partner_res not in self.seen_res:
                            self.seen_res[r] = 1
                            ss = "("
                        elif partner_res in self.seen_res:
                            if self.seen_res[partner_res] > 1:
                                ss = "."
                            else:
                                ss = ")"
                                self.seen_res[r] = 1
                                self.seen_res[partner_res] += 1
                                break
                    elif r not in self.seen_res:
                        ss = "."

                if saved_bp is not None:
                    self.seen_bp[saved_bp] = 1

                ss_res.append(secondary_structure.Residue(r.name, ss, r.num,
                                                          r.chain_id, r.uuid, r.i_code))
            ss_chains.append(secondary_structure.Chain(ss_res))
            best_chain = self._get_next_chain(motif)

            if best_chain is None:
                break
            self.chains.remove(best_chain)
            self.open_chains.append(best_chain)

        ss = secondary_structure.SecondaryStructure(chains=ss_chains)
        self._setup_basepairs_and_ends(ss, motif)

        return ss


class StructureSecondaryFactory(object):
    def __init__(self):
        self.parser = MotiftoSecondaryStructure()

    def _get_basepairs(self, sstree, ss):
        basepairs, ends = [], []
        for n in sstree:
            if n.data.type == ss_tree.SS_Type.SS_BP:
                res1 = n.data.ss_chains[0].residues[0]
                res2 = n.data.ss_chains[1].residues[0]

                bp = secondary_structure.Basepair(res1, res2)
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
                if n.parent.data.type == ss_tree.SS_Type.SS_SEQ_BREAK:
                    ends.append(bp)

        return basepairs, ends

    def _get_elements(self, sstree, ss):
        elements = { 'ALL' : [] }

        for n in sstree:
            if n.parent_index() == -1 or n.parent.data.type == ss_tree.SS_Type.SS_SEQ_BREAK:
                continue

            if n.data.type == ss_tree.SS_Type.SS_BP and \
               n.parent.data.type == ss_tree.SS_Type.SS_BP:
                nodes = [n,n.parent]
                bp_nodes = [n,n.parent]
                type_name = "BP_STEP"

            elif n.data.type != ss_tree.SS_Type.SS_SEQ_BREAK and \
                n.parent.data.type == ss_tree.SS_Type.SS_BP:
                nodes = [n.parent, n]
                bp_nodes = [n.parent]
                type_name = n.data.what()[3:]
                for c in n.children:
                    if c is None:
                        continue
                    if c.data.type == ss_tree.SS_Type.SS_BP or \
                       c.data.type == ss_tree.SS_Type.SS_PSEUDO_BP:
                        nodes.append(c)
                        bp_nodes.append(c)
                    else:
                        print n.data.what(), c.data.what()
                        raise ValueError("unexpected connectivity in ss_tree")
            else:
                continue


            chains = self._get_chains(ss, nodes)
            ends = []
            for bp_n in bp_nodes:
                bp = ss.get_bp(bp_n.data.ss_chains[0].residues[0],
                               bp_n.data.ss_chains[1].residues[0])
                ends.append(bp)
            e = secondary_structure.SecondaryStructureMotif(type_name, ends, chains)
            bps = []
            res = {r : 1 for r in e.residues() }
            for bp in ss.basepairs:
                if bp.res1 in res and bp.res2 in res:
                    bps.append(bp)
            e.basepairs = bps
            if type_name not in elements:
                elements[type_name] = []
            end_ids = []
            for end in e.ends:
                end_ids.append(secondary_structure.assign_end_id(e, end))
            e.end_ids = end_ids
            elements[type_name].append(e)
            elements['ALL'].append(e)
        return elements

    def _get_chains(self, ss, nodes):
        res = []
        for n in nodes:
            for c in n.data.ss_chains:
                res.extend(c.residues)
        res.sort(key=lambda x: x.num)
        chains = []
        c_res = []
        last = -1
        for r in res:
            is_chain_start = 0
            for c in ss.chains:
                if c.first() == r:
                    is_chain_start = 1
                    break

            if last == -1:
                pass
            elif last+1 != r.num or is_chain_start:
                chains.append(secondary_structure.Chain(c_res))
                c_res = []
            c_res.append(r)
            last = r.num
        if len(c_res) > 0:
            chains.append(secondary_structure.Chain(c_res))
        return chains

    def get_structure(self, sequence=None, dot_bracket=None, base_ss=None):
        if   sequence is not None and dot_bracket is not None:
            sstree = ss_tree.SS_Tree(sequence, dot_bracket)
        elif base_ss is not None:
            sstree = ss_tree.SS_Tree(ss=base_ss)
        else:
            raise ValueError("supply sequence and dot_bracket strings or a" + \
                             " SecondaryStructure object")
        ss = sstree.ss

        ss.basepairs, ss.ends = self._get_basepairs(sstree, ss)
        ss.elements = self._get_elements(sstree, ss)

        return ss

    def secondary_structure_from_motif(self, m):
        ss = self.parser.to_secondary_structure(m)
        self.parser.reset()
        return ss


def ss_id_to_seq_and_db(ss_id):
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
    return seq, ss


def ss_id_to_ss_tree(ss_id):
    seq, ss = ss_id_to_seq_and_db(ss_id)
    return ss_tree.SS_Tree(ss, seq)


def ss_id_to_secondary_structure(ss_id):
    seq, ss = ss_id_to_seq_and_db(ss_id)
    return factory.get_structure(seq, ss)


factory = StructureSecondaryFactory()
