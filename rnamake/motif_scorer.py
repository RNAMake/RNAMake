import util

class MotifScorer(object):

    def __init__(self):
        self.bp_ref_energy = self._bp_reference_energy_table()
        self.unpaired_pentalty = 4.0

    def _score_cWW_bp(self, bp):
        if not ( util.wc_bp(bp) or util.gu_bp(bp) ):
            return 2
        bpstr = bp.res1.rtype.name[0]+bp.res2.rtype.name[0]
        if   bpstr == "GC" or bpstr == "CG":
            return -2
        elif bpstr == "AU" or bpstr == "UA":
            return -1
        else:
            return -0.5

    def score(self, m):
        score = 0
        for bp in m.basepairs:
            if   bp.bp_type == "cW-W":
                score += self._score_cWW_bp(bp)
            elif bp.bp_type in self.bp_ref_energy:
                score += self.bp_ref_energy[bp.bp_type]
            else:
                score += 6.1122

        for r in m.residues():
            bps = m.get_basepair(uuid1=r.uuid)
            if len(bps) == 0:
                score += self.unpaired_pentalty

        score += len(m.residues())*.05
        return score

    def _score_cWW_bp_new(self, bp, s):
        if not ( util.wc_bp(bp, s) or util.gu_bp(bp, s)):
            return 2
        res1 = s.get_residue(uuid=bp.res1_uuid)
        res2 = s.get_residue(uuid=bp.res2_uuid)
        bpstr = res1.short_name()+res2.short_name()
        if   bpstr == "GC" or bpstr == "CG":
            return -2
        elif bpstr == "AU" or bpstr == "UA":
            return -1
        else:
            return -0.5

    def score_elements(self, s, basepairs):
        score = 0
        for bp in basepairs:
            if   bp.bp_type == "cW-W":
                score += self._score_cWW_bp_new(bp, s)
            elif bp.bp_type in self.bp_ref_energy:
                score += self.bp_ref_energy[bp.bp_type]
            else:
                score += 6.1122

        for r in s:
            found = 0
            for bp in basepairs:
                if r.uuid == bp.res1_uuid or r.uuid == bp.res2_uuid:
                    found = 1
                    break
            if not found:
                score += self.unpaired_pentalty

        score += s.num_residues()*.05
        return score

    def _bp_reference_energy_table(self):
        bp_ref_energy = {
            'cm-':6.1122307919,
            'cM-M':6.1122307919,
            'tW+W':3.11366762268,
            'c.+M':5.68522906698,
            '.W+W':6.1122307919,
            'tW-M':2.42283130036,
            'tm-M':2.71577524698,
            'cW+M':3.33339125508,
            '.W-W':4.33166562348,
            'cM+.':6.1122307919,
            'c.-m':6.1122307919,
            'cM+W':4.4042238922,
            'tM+m':6.1122307919,
            'tM-W':3.02141948251,
            'cm-m':5.12076349023,
            'cM-W':6.1122307919,
            'cW-W':0.056986982519,
            'c.-M':5.43544907015,
            'cm+M':2.7132962365,
            'cm-M':3.23361276018,
            '....':4.18066203386,
            'cm-W':4.36687710812,
            'tM-m':2.83911913314,
            'c.-W':6.1122307919,
            'cM+m':5.68522906698,
            'cM-m':3.12321871743
        }

        return bp_ref_energy

