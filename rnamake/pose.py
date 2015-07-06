import motif
import motif_type
import x3dna
import eternabot.sequence_designer as sequence_designer
import secondary_structure

class Pose(motif.Motif):
    """
    Poses are for loading for RNA structures that have multiple motifs.
    Like motifs poses store 3D coordinates as well as 3DNA parameters but
    in addition also contains information about what types of motifs are
    contained in the structure as well as the ability to generate newly
    designed sequences.

    :param mdir: the path to a pose directory that contains required files
    :type mdir: str

    :param pdb: the path to a pdb file to create this pose object, will also
        creates ref_frames.dat file and dssr out file in current directory
    :type pdb: str

    Attributes
    ----------
    `base_pairs` : List of Basepair objects
        All the basepair info determined from 3DNA
    `beads` : List of Bead objects
        All the beads in the 3 bead residue model for all the residues in
        structure object
    `dir` : str
        full path to directory
    `ends` : List of the Basepair objects
        that are at the end of Pose together its not necessarily the
        first and last Basepairs
    `structure` : Structure object
        holds 3D coordinate data
    `name` : str
        the name of the directory without entire path
    `motifs` : Motifs object
        contains all the motifs found in this RNA Pose
    `designable` : Dict object
        stores whether a given basepair can be optimized in design or not

    """
    def __init__(self, mdir=None, pdb=None):
        self.structure = None
        self.mdir, self.name, self.ends = "", "", []
        self.beads, self.score, self.basepairs = [], 0, []
        self.designable = {}
        self.mtype = motif_type.UNKNOWN
        self.end_ids = []
        self.secondary_structure = secondary_structure.SecondaryStructure()
        #self._setup(mdir, pdb)
        if pdb is not None or mdir is not None:
            self._setup_motifs()

    def _setup_motifs(self):
        """
        parses x3dna output to find all motifs and characterizes them in a
        Motifs object, also goes through all helical motifs and checks to see
        if their basepairs can be designed in sequence optimization
        """
        x = x3dna.X3dna()
        motifs = x.get_motifs(self.mdir+"/"+self.name)
        self.motifs = Motifs(motifs, self)
        basepairs = []
        for helix in self.helices():
            basepairs.extend(helix.basepairs)
        for bp in basepairs:
            found_bp = 0
            for m in self.all_motifs():
                if m.mtype == motif_type.HELIX:
                    continue
                found = m.get_basepair(bp_uuid=bp.uuid)
                if len(found) > 0:
                    found_bp = 1
            if not found_bp:
                self.designable[bp.uuid]=1

    def twoways(self):
        """
        wrapper for motifs.twoway, returns all twoway junctions motif object
        found in this pose
        """
        return self.motifs.twoways

    def nways(self):
        """
        wrapper for motifs.nway, returns all nway junctions motif object
        found in this pose
        """
        return self.motifs.nways

    def hairpins(self):
        """
        wrapper for motifs.hairpins, returns all hairpin motif object
        found in this pose
        """
        return self.motifs.hairpins

    def helices(self):
        """
        wrapper for motifs.helices, returns all helix motif objects
        found in this pose
        """
        return self.motifs.helices

    def single_strands(self):
        """
        wrapper for motifs.single_strands, returns all single stranded motif object
        found in this pose
        """
        return self.motifs.single_strands

    def all_motifs(self):
        """
        wrapper for motifs.all_motifs, returns all motif objects
        found in this pose
        """
        return self.motifs.all_motifs

    def designable_sequence(self):
        """
        returns sequence containing Ns for locations that can be changed,
        these spots correspond with helix areas that most likely will not
        affect the overall fold
        """
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

    def optimized_sequence(self, ss=None):
        """
        returns the eternabot optimized sequence for the pose
        """
        if len(self.chains()) > 2:
            raise ValueError("cannot get optimized sequence with more then 2 "\
                             "chains cannot call RNAFold or RNAcoFold")

        seq = self.designable_sequence()
        if ss == None:
            ss  = self.dot_bracket()

        designer = sequence_designer.SequenceDesigner()
        results = designer.design(ss, seq)
        score = results[0].score
        return results[0].sequence



#TODO implement finding tertiary contacts
class Motifs(object):
    def __init__(self, x3dna_motifs=None, p=None):
        self.twoways, self.nways, self.hairpins, self.helices = [], [], [], []
        self.single_strands, self.all_motifs = [], []
        if x3dna_motifs is not None:
            self._setup_from_x3dna_motif(x3dna_motifs, p)

    def _setup_from_x3dna_motif(self, x3dna_motifs, p):
        for xm in x3dna_motifs:
            m = self._convert_x3dna_to_motif(xm, p)
            self.all_motifs.append(m)
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
        m.setup_basepair_ends()
        return m


