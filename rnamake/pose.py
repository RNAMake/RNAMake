import motif
import eternabot.sequence_designer as sequence_designer


class Pose(motif.Motif):
    #TODO finish implementing
    def __init__(self, mdir=None, pdb=None):
        self.structure = None
        self.mdir, self.name, self.ends = "", "", []
        self.beads, self.score, self.basepairs = [], 0, []
        self.designable = {}
        self._setup(mdir, pdb)

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




