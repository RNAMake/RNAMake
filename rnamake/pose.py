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
    def __init__(self):
        super(self.__class__, self).__init__()
        self.designable = {}
        self.motif_dict = {
            motif_type.TWOWAY        : [],
            motif_type.NWAY          : [],
            motif_type.HAIRPIN       : [],
            motif_type.HELIX         : [],
            motif_type.SSTRAND       : [],
            motif_type.TCONTACT      : [],
            motif_type.ALL           : []
        }

    def __repr__(self):
        """
        is called when motif is printed
        """
        return "<Pose(\n\tstructure='%s', \n\tends='%s')>" % (
        self.structure,len(self.ends))

    def motifs(self, mtype):
        return self.motif_dict[mtype]

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
        return results[0].sequence

