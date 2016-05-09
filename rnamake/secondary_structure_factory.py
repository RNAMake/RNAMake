import secondary_structure
import secondary_structure_parser
import util


class MotiftoSecondaryStructure(object):
    def __init__(self):
        self.chains = []
        self.open_chains = []
        self.seen_res = {}
        self.seen_bp = {}

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

    def _setup_basepairs_and_ends(self, struct, motif):
        ss_bps = []
        for bp in self.seen_bp.keys():
            res1 = struct.get_residue(uuid=bp.res1.uuid)
            res2 = struct.get_residue(uuid=bp.res2.uuid)
            ss_bps.append(secondary_structure.Basepair(res1, res2, bp.uuid))
        ss_ends = []
        for end in motif.ends:
            res1 = struct.get_residue(uuid=end.res1.uuid)
            res2 = struct.get_residue(uuid=end.res2.uuid)
            end_bp = None
            for bp in ss_bps:
                if bp.res1 == res1 and bp.res2 == res2:
                    end_bp = bp
                    break
            if end_bp is None:
                motif.to_pdb("test.pdb")
                raise ValueError("did not properly find end in generating ss")
            ss_ends.append(end_bp)

        return secondary_structure.Motif(struct, ss_bps, ss_ends)

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

                ss_res.append(secondary_structure.Residue(r.short_name(), ss, r.num,
                                                          r.chain_id, r.uuid, r.i_code))
            ss_chains.append(secondary_structure.Chain(ss_res))
            best_chain = self._get_next_chain(motif)

            if best_chain is None:
                break
            self.chains.remove(best_chain)
            self.open_chains.append(best_chain)

        struct = secondary_structure.Structure(chains=ss_chains)
        m = self._setup_basepairs_and_ends(struct, motif)

        return m


class StructureSecondaryFactory(object):
    def __init__(self):
        self.parser = MotiftoSecondaryStructure()

    def secondary_structure_from_motif(self, m):
        ss = self.parser.to_secondary_structure(m)
        self.parser.reset()
        return ss

    def get_helices(self, ss):
        bp_steps = ss.motifs('BP_STEP')
        seen = {}
        groups = []
        seen [ bp_steps[0] ] = 1
        group = [ bp_steps[0] ]

        while len(seen) < len(bp_steps):
            if len(group) == 0:
                for bp_step in bp_steps:
                    if bp_step not in seen:
                        group.append(bp_step)
                        break

            found = 0
            ends = []
            for step in group:
                ends.extend(step.ends)

            for bp_step in bp_steps:
                if bp_step in seen:
                    continue

                for end in ends:
                    if end in bp_step.ends:
                        found = 1
                        break

                if found:
                    seen[bp_step] = 1
                    group.append(bp_step)
                    break

            if not found:
                groups.append(group)
                group = []

        if len(group) > 0:
            groups.append(group)


        helices = []
        for g in groups:
            chains = self._get_ss_chains(ss, g)
            e = secondary_structure.SecondaryStructureMotif("HELIX", [], chains)
            helices.append(e)

        return helices

    def motif(self, sequence=None, dot_bracket=None):
        parser = secondary_structure_parser.SecondaryStructureParser()
        return parser.parse_to_motif(sequence, dot_bracket)

    def pose(self, sequence=None, dot_bracket=None):
        parser = secondary_structure_parser.SecondaryStructureParser()
        return parser.parse_to_pose(sequence, dot_bracket)


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


factory = StructureSecondaryFactory()
