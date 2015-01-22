import motif
import motif_type
import x3dna
import eternabot.sequence_designer as sequence_designer


class Pose(motif.Motif):
    #TODO finish implementing
    def __init__(self, mdir=None, pdb=None):
        self.structure = None
        self.mdir, self.name, self.ends = "", "", []
        self.beads, self.score, self.basepairs = [], 0, []
        self.designable = {}
        self._setup(mdir, pdb)
        self._setup_motifs()

    def _setup_motifs(self):
        x = x3dna.X3dna()
        motifs = x.get_motifs(self.mdir+"/"+self.name)
        self.motifs = Motifs(motifs, self)

    def twoways(self):
        return self.motifs.twoways

    def nways(self):
        return self.motifs.nways

    def hairpins(self):
        return self.motifs.hairpins

    def helices(self):
        return self.motifs.helices

    def single_strands(self):
        return self.motifs.single_strands

    def designable_sequence(self):
        seq = ""
        for c in self.chains():
            for r in c.residues:
                bps = self.get_basepair(res1=r)
                s = r.rtype.name[0]
                for bp in bps:
                    if bp.uuid in self.designable:
                        s = "N"
                        break
                seq += s
            seq += "&"
        return seq[:-1]

    def optimized_sequence(self):
        if len(self.chains()) > 2:
            raise ValueError("cannot get optimized sequence with more then 2 "\
                             "chains cannot call RNAFold or RNAcoFold")

        seq = self.designable_sequence()
        ss  = self.secondary_structure()

        designer = sequence_designer.SequenceDesigner()
        results = designer.design(ss, seq)
        return results[0]['end'][0]

class Motifs(object):
    def __init__(self, x3dna_motifs, p):
        self.twoways, self.nways, self.hairpins, self.helices = [], [], [], []
        self.single_strands = []
        for xm in x3dna_motifs:
            m = self._convert_x3dna_to_motif(xm, p)
            self._assign_motif_by_type(m)

    def _assign_motif_by_type(self, m):
        if   m.mtype == motif_type.TWOWAY:
            self.twoways.append(m)
        elif m.mtype == motif_type.NWAY:
            self.nways.append(m)
        elif m.mtype == motif_type.HAIRPIN:
            self.hairpins.append(m)
        elif m.mtype == motif_type.HELIX:
            self.helices.append(m)
        else:
            self.single_strands.append(m)

    def _convert_x3dna_to_motif(self, xm, p):
        res = []
        for xr in xm.residues:
            r = p.get_residue(num=xr.num, chain_id=xr.chain_id, i_code=xr.i_code)
            res.append(r)
        basepairs = []
        for r in res:
            bps = p.get_basepair(res1=r)
            for bp in bps:
                if bp.res1 in res and bp.res2 in res and bp not in basepairs:
                    basepairs.append(bp)
        m = motif.Motif(mtype=xm.mtype)
        m.structure._build_chains(res)
        m.structure._cache_coords()
        m.basepairs = basepairs
        m._cache_basepair_frames()
        return m


