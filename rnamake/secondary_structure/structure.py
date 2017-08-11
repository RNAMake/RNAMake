import uuid

# rnamake imports
from rnamake import motif_type
from rnamake import exceptions
from rnamake import primitives
from rnamake import motif_type_directed_graph

class Residue(primitives.Residue):
    """
    An extremely stripped down container object for use for keeping track of
    secondary structure using dot bracket notation. Dot bracket notation is
    represented by a string of equal length of a sequence, example:

    sequence:    'GCAAAACG'\n
    dot bracket: '((....))'

    '(' represent the 5' end of a basepair, ')' the 3' end of the same base
    pair. Thus there should be an equal number of '(' and ')' symbols. '.'
    represents an unpaired residue.

    This class along with others in this module are generally
    not instatinated outside RNAStructure,  Motif  and Pose.

    :param name: name of residue, etc A, G, C, U
    :param dot_bracket: dot bracket notation for secondary structure either '(', '.', ')'
    :param num: number of residue
    :param chain_id: chain id of residue, etc "A" or "B"
    :param uuid: residue unique indentifier
    :param i_code: residue insertaion code, usually ""

    :type name: str
    :type dot_bracket: str
    :type num: int
    :type chain_id: str
    :type uuid: uuid.uuid1
    :type i_code: str

    :attributes:

    `name`: str
        name of residue, etc A, G, C, U
    `dot_bracket`: str
        dot bracket notation for secondary structure either '(', '.', ')'
    `num`: int
        number of residue
    `chain_id`: str
        chain id of residue, etc "A" or "B"
    `uuid`: uuid.uuid1
        residue unique indentifier
    `i_code`: str
        residue insertaion code, usually ""

    :examples:

    ..  code-block:: python

        # create a new residue
        >>> r = Residue("G", "(", 10, "A", uuid.uuid1())

    """

    __slots__= [
        "_name",
        "_dot_bracket",
        "_num",
        "_chain_id",
        "_uuid",
        "_i_code"]

    def __init__(self, name, dot_bracket, num, chain_id, i_code=None, r_uuid=None):
        self._dot_bracket = dot_bracket
        super(self.__class__, self).__init__(name, num, chain_id, i_code, r_uuid)

    @classmethod
    def from_str(cls, s):
        spl = s.split(",")
        return cls(spl[0], spl[1], int(spl[2]), spl[3], spl[4])

    @classmethod
    def copy(cls, r, new_uuid=0):
        """
        creates copy of current residue

        :return: copy of instatnce
        :rtype: secondary_structure.Residue
        """
        r_uuid = r._uuid
        if new_uuid:
            r_uuid = uuid.uuid1()

        return cls(r._name, r._dot_bracket, r._num,
                   r._chain_id, r._i_code, r_uuid)

    def __repr__(self):
        return "<SecondaryStructureResidue('%s%d%s chain %s')>" % (
            self._name, self._num, self._i_code, self._chain_id)

    def get_str(self):
        """
        stringify residue object. can be converted back with
        :func:`str_to_residue`

        :return: stringifed verision of residue
        :rtype: str
        """

        return self._name + "," + self._dot_bracket + "," + str(self._num) + "," + \
               str(self._chain_id) + "," + str(self._i_code)

    @property
    def dot_bracket(self):
        return self._dot_bracket

    def set_name(self, n):
        self._name = n

    def short_name(self):
        return self._name


class Chain(primitives.Chain):
    """
    secondary structure chain container. Contains a chain of connected
    residues. Chain should be from 5' to 3'.

    :param residues: the residues that are to be included in chain. Optional.
    :type residues: list of Residues.

    :attributes:
    `residues` : list of Residues
        Residues contained in current chain

    """

    __slots__ = ["_residues"]

    def __init__(self, residues):
        super(self.__class__, self).__init__(residues)

    @classmethod
    def from_str(cls, s):
        """
        converts a chain from string generated from :func:`Chain.to_str`

        :param s: string created by Chain.to_str()
        :type s: str

        :return: chain from str
        :rtype: secondary_structure.Chain
        """

        spl = s.split(";")
        residues = []
        for r_str in spl[:-1]:
            r = Residue.from_str(r_str)
            residues.append(r)
        return cls(residues)

    @classmethod
    def copy(cls, c):
        """
        creates deep copy of chain instance

        :return: copy of chain
        :rtype: secondary_structure.Chain

        """

        residues = []
        for r in c:
            residues.append(Residue.copy(r))
        return cls(residues)

    def __repr__(self):
        seq = ""
        for r in self._residues:
            seq += r.name

        return "<SecondaryStructureChain( " + seq + ")"

    def sequence(self):
        """
        gets the string verision of the sequence in this chain

        :returns: string of sequence of chain
        :rtype: str

        :examples:

        ..  code-block:: python

            >>> from rnamake.unittests import instances
            >>> c = instances.secondary_structure_chain()

            # see the residue of each residue in the chain, all Gs
            >>> for r in c: print r
            ...
            <SecondaryStructureResidue('G13 chain A')>
            <SecondaryStructureResidue('G14 chain A')>
            <SecondaryStructureResidue('G15 chain A')>
            <SecondaryStructureResidue('G16 chain A')>
            <SecondaryStructureResidue('G17 chain A')>
            <SecondaryStructureResidue('G18 chain A')>
            <SecondaryStructureResidue('G19 chain A')>
            <SecondaryStructureResidue('G20 chain A')>
            <SecondaryStructureResidue('G21 chain A')>
            <SecondaryStructureResidue('G22 chain A')>
            <SecondaryStructureResidue('G23 chain A')>
            <SecondaryStructureResidue('G24 chain A')>

            >>> c.sequence()
            u'GGGGGGGGGGGG'
        """

        seq = ""
        for r in self._residues:
            seq += r.name
        return seq

    def dot_bracket(self):
        """
        gets the string verision of the secondary structure in
        dot_bracket notation of this chain

        :returns: string of secondary structure in dot bracket notation of
            chain
        :rtype: str

        :examples:

        ..  code-block:: python

            >>> from rnamake.unittests import instances
            >>> c = instances.secondary_structure_chain()
            >>> c.dot_bracket()
            u'(((((((((((('
        """

        db = ""
        for r in self._residues:
            db += r.dot_bracket
        return db

    def get_str(self):
        """
        stringify chain object. can be converted back with
        :func:`str_to_chain`

        :return: stringifed verision of chain
        :rtype: str
        """

        s = ""
        for r in self._residues:
            s += r.to_str() + ";"
        return s


class Structure(primitives.Structure):
    """
    lightweight container class for storing secondary structure information
    for an entire RNA.

    :param chains: secondary_structure.Chains that belong to this structure,
        this is done when a structure is being built from an existing 3D
        structure.Structure instance
    :param sequence: sequence of RNA of interest, e.g. "AAAGGGCCC",
    :param dot_bracket: dot bracket notation of the secondary structure of
        RNA of interes, e.g. "(((())))"

    :type chains: list of secondary_structure.Chains
    :type sequence: str
    :type dot_bracket: str

    :attributes:

    `chains`: list of secondary_structure.Chains
        the chains of RNA residues in this structure

    :examples:

    ..  code-block:: python

        # create a new structure
        >>> from rnamake import secondary_structure
        >>> s = secondary_structure.Structure(sequence="GCGAAAACGC",
                                              dot_bracket="(((....)))")

        >>> print s.sequence()
        GCGAAAACGC
        >>> print s.dot_bracket()
        (((....)))

        >>> s.get_residue(num=1)
        <SecondaryStructureResidue('G1 chain A')>
    """

    __slots__ = [
        "_chain_cuts",
        "_residues"
    ]

    def __init__(self, residues, chain_cuts):
        super(self.__class__, self).__init__(residues, chain_cuts)

    @classmethod
    def from_str(cls, s):
        """
        converts a structure from string generated from :func:`Structure.to_str`

        :param s: string created by Structure.to_str()
        :type s: str

        :return: structure from str
        :rtype: secondary_structure.Structure
        """
        spl = s.split("|")
        chains = []
        for c_str in spl[:-1]:
            c = Chain.from_str(c_str)
            chains.append(c)
        return cls(chains)

    @classmethod
    def copy(cls, s):
        """
        creates a deep copy of structure instance

        :return: copy of structure
        :rtype: secondary_structure.Structure
        """
        new_chains = [ Chain.copy(c) for c in s]
        return cls(new_chains)

    def sequence(self):
        """
        Concats the sequence of each Chain into one sequence for the entire
        RNA

        :return: sequence of structure
        :rtype: seq
        """
        sequences = [x.sequence() for x in self.get_chains()]
        return "&".join(sequences)

    def dot_bracket(self):
        """
        Concats the secondary structure in the form of dot bracket notation
        of each Chain into one sequence for the entire RNA

        :return: sequence of structure
        :rtype: seq
        """
        dot_brackets = [x.dot_bracket() for x in self.get_chains()]
        return "&".join(dot_brackets)

    def get_str(self):
        """
        generates a stringified verision of this instance.

        :returns: stringified verision of structure
        :rtype: str
        """
        s = ""
        for c in self._chains:
            s += c.to_str() + "|"
        return s

    def get_chains(self):
        pos = 0
        res = []
        chains = []
        for i, r in enumerate(self._residues):
            if self._chain_cuts[pos] == i:
                c = Chain(res)
                chains.append(c)
                res = [r]
                pos += 1
            else:
                res.append(r)

        if len(res) > 0:
            chains.append(Chain(res))
        return chains


class Basepair(primitives.Basepair):
    """
    :param res1: First residue in basepair
    :param res2: Second residue in basepair
    :param bp_uuid: basepair unique indentifier

    :type res1: secondary_structure.Residue
    :type res2: secondary_structure.Residue
    :type bp_uuid: uuid.uuid1

    :attributes:

    `res1` : secondary_structure.Residue
        First residue in basepair
    `res2` : secondary_structure.Residue
        Second residue in basepair
    `uuid`: uuid.uuid1
        unique id to indentify this basepair when locating it in a motif or
        pose
    """

    __slots__ = [
        "_res1_uuid",
        "_res2_uuid",
        "_name",
        "_uuid",
        "_bp_type"
    ]

    def __init__(self, res1_uuid, res2_uuid, name, bp_uuid=None):
        self._res1_uuid = res1_uuid
        self._res2_uuid = res2_uuid
        self._name = name
        self._uuid = bp_uuid
        self._bp_type = "cW-W"
        if self._uuid is None:
            self._uuid = uuid.uuid1()

    @classmethod
    def copy(cls, bp):
        return cls(bp._res1_uuid, bp._res2_uuid, bp._name, bp._uuid)

    def __repr__(self):
        return "<SecondaryStructureBasepair("+self._name+")>"

    def partner(self, r_uuid):
        """
        get the other basepairing partner of a residue will throw an error
        if the supplied residue is not contained within this basepair

        :param res: the residue that you want to get the partner of
        :type res: secondary_structure.Residue object
        """

        if   r_uuid == self._res1_uuid:
            return self.res2_uuid
        elif r_uuid == self._res2_uuid:
            return self.res1_uuid
        else:
            raise exceptions.SecondaryStructureException(
                "call partner with a residue not in basepair")

    @property
    def name(self):
        return _name

    @property
    def res1_uuid(self):
        return self._res1_uuid

    @property
    def res2_uuid(self):
        return self._res2_uuid

    @property
    def uuid(self):
        return self._uuid

    @property
    def bp_type(self):
        return self._bp_type


class RNAStructure(primitives.RNAStructure):
    """
    Complete secondary structure container for representing a RNA. Contains
    both the sequence indentity of each residue with its corresponding dot
    bracket notation symbol but also includes basepair objects to represent the
    pairs between residues. This class parallels the
    rna_structure.RNAStructure class for describing RNA with 3D coordinates.
    Having a parallel class for secondary structure makes it simple to move
    between secondary structure and full atom representations of RNA.
    RNAStructure is rarely called directly but serves as an abstract class
    for both Motif and Pose so for example use please see those classes.

    :param structure: structure containing residue and chain information for
        this RNAStructure instance
    :param basepairs: basepairs contained in RNAStructure
    :param ends: the basepairs at the end of chains. These define connection
        points to other RNAStructures
    :param name: name of RNAStructure
    :param path: location of where RNAStructure originated from, this is just
        a place holder for converting from rna_structure.RNAStructure
    :param score: the score generated by motif_scorer.MotifScorer
    :param end_ids: strings indenifying the secondary structure and sequence
        in the perspective of a given basepair end see
        :func:`assign_end_id_new`

    :type structure: secondary_structure.Structure
    :type basepairs: list of secondary_structure.Basepairs
    :type ends: list of secondary_structure.Basepairs
    :type name: str
    :type path: str
    :type score: float
    :type end_ids: list of strs

    :attributes:

    `structure`: secondary_structure.Structure
        structure containing residue and chain information for
        this RNAStructure instance
    `basepairs`: list of secondary_structure.Basepairs
        Basepairs between residues
    `ends`: list of secondary_structure.Basepairs
        Basepair ends where RNA structures can be connected
    `name`: str
        the name of the RNAStructure
    `path`: str
        location of where RNAStructure originated from, this is just
        a place holder for converting from rna_structure.RNAStructure
    `score` : float
        the score generated by motif_scorer.MotifScorer, estimates secondary
        structure stability
    `end_ids`: list of strs
        strings indenifying the secondary structure and sequence
        in the perspective of a given basepair end see
        :func:`assign_end_id_new`

    """

    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_end_ids"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, name):
       super(self.__class__, self).__init__(structure, basepairs, ends,
                                            end_ids, name)

    @classmethod
    def from_str(cls, s):
        pass

    @classmethod
    def copy(cls, rs):
        """
        creates deep copy of RNAStructure instance

        :returns: copy of instance
        :rtype: RNAStructure
        """

        n_ss = Structure.copy(rs._structure)
        basepairs, ends = [], []
        for bp in rs._basepairs:
            basepairs.append(Basepair.copy(bp))
        for end in rs._ends:
            i = rs._basepairs.index(end)
            ends.append(basepairs[i])

        return cls(n_ss, basepairs, ends, rs._end_ids[::])

    def __repr__(self):
        return "<secondary_structure.RNAStructure( " + self.sequence() + " " +\
                    self.dot_bracket() + " )"

    def sequence(self):
        """
        wrapper for :func:`Structure.sequence`
        """

        return self._structure.sequence()

    def dot_bracket(self):
        """
        wrapper for :func:`Structure.dot_bracket`
        """

        return self._structure.dot_bracket()

    def replace_sequence(self, seq):
        """
        changes the sequence of structure.

        :param seq: the new sequence of structure
        :type seq: str

        :returns: None
        """

        spl = seq.split("&")
        seq2 = "".join(spl)

        if len(seq2) != self._structure.num_residues():
            raise exceptions.SecondaryStructureException(
                "cannot replace sequence, sequence length is different then "
                "the number of residues")

        for i, r in enumerate(self._structure.iter_res()):
            r.set_name(seq2[i])

    def get_str(self):
        pass

    def update(self):
        for i, end in enumerate(self._ends):
            self._end_ids[i] = primitives.rna_structure.assign_end_id(
                                    self._structure, self._basepairs, self._ends, end)


class Motif(RNAStructure):
    """
    Complete secondary structure container for representing an RNA Motif.
    Contains both the sequence indentity of each residue with its corresponding
    dot bracket notation symbol but also includes basepair objects to represent
    the pairs between residues. This class parallels the motif.Motif class for
    describing RNA with 3D coordinates.
    Having a parallel class for secondary structure makes it simple to move
    between secondary structure and full atom representations of RNA.

    :param structure: structure containing residue and chain information for
        this RNAStructure instance
    :param basepairs: basepairs contained in RNAStructure
    :param ends: the basepairs at the end of chains. These define connection
        points to other Motifs
    :param name: name of Motif
    :param path: location of where Motif originated from, this is just
        a place holder for converting from motif.Motif
    :param mtype: motif_type enum value, to indentify what type of motif this
        is
    :param score: the score generated by motif_scorer.MotifScorer
    :param end_ids: strings indenifying the secondary structure and sequence
        in the perspective of a given basepair end see
        :func:`assign_end_id_new`
    :param id: unique indentifer for motif
    :param r_struct: an RNAStructure instance to setup from

    :type structure: secondary_structure.Structure
    :type basepairs: list of secondary_structure.Basepairs
    :type ends: list of secondary_structure.Basepairs
    :type name: str
    :type path: str
    :type mtype: motif_type
    :type score: float
    :type end_ids: list of strs
    :type id: uuid.uuid1
    :type r_structure: secondary_structure.RNAStructure

    :attributes:

    `structure`: secondary_structure.Structure
        structure containing residue and chain information for
        this Motif instance
    `basepairs`: list of secondary_structure.Basepairs
        Basepairs between residues
    `ends`: list of secondary_structure.Basepairs
        Basepair ends where RNA structures can be connected
    `name`: str
        the name of the Motif
    `path`: str
        location of where Motif originated from, this is just
        a place holder for converting from motif.Motif
    `score` : float
        the score generated by motif_scorer.MotifScorer, estimates secondary
        structure stability
    `end_ids`: list of strs
        strings indenifying the secondary structure and sequence
        in the perspective of a given basepair end see
        :func:`assign_end_id_new`

    """

    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_end_ids",
        "_mtype",
        "_uuid",
        "_name",
        "_block_end_add"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, mtype, name,
                 m_uuid=None):
        super(RNAStructure, self).__init__(structure, basepairs,
                                           ends, end_ids, name)
        self._mtype = mtype
        self._uuid = m_uuid
        if self._uuid is None:
            self._uuid = uuid.uuid1()

    def __repr__(self):
        return "<secondary_structure.Motif( " + self.sequence() + " " + self.dot_bracket() + " )"

    def copy(self):
        """
        creates a deep copy of Motif instance

        :returns: deep copy of instance
        :rtype: secondary_structure.Motif
        """

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
        """
        creates a deep copy of Motif instance replacing residues and basepairs
        that were already created. This is used when a pose is being copied to
        make sure there are duplicate copies of residues and basepairs.

        :param res: list of residues that were already copied
        :param bps: list of basepairs that were already copied

        :type res: list of secondary_structure.Residues
        :type bps: list of secondary_structure.Basepairs

        :returns: deep copy of instance
        :rtype: secondary_structure.Motif
        """

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
        """
        creates a stringified verision of this instance

        :returns: stringified verision of instance
        :rtype: str
        """

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

    @property
    def uuid(self):
        return self._uuid


class Aligner(primitives.Aligner):
    def __init__(self):
        super(self.__class__, self).__init__()

    def get_aligned_motif(self, ref_bp, m):
        return m

    def align(self, ref_bp, m):
        return m


class Pose(RNAStructure):
    """
    Complete secondary structure container for representing an RNA Pose.
    Contains both the sequence indentity of each residue with its corresponding
    dot bracket notation symbol but also includes basepair objects to represent
    the pairs between residues. This class parallels the pose.Pose class for
    describing RNA with 3D coordinates.
    Having a parallel class for secondary structure makes it simple to move
    between secondary structure and full atom representations of RNA.
    A pose is a composite structure containing more then one motif.

    :param structure: structure containing residue and chain information for
        this RNAStructure instance
    :param basepairs: basepairs contained in RNAStructure
    :param ends: the basepairs at the end of chains. These define connection
        points to other Motifs
    :param name: name of Motif
    :param path: location of where Pose originated from, this is just
        a place holder for converting from pose.Pose

    :param score: the score generated by motif_scorer.MotifScorer
    :param end_ids: strings indenifying the secondary structure and sequence
        in the perspective of a given basepair end see
        :func:`assign_end_id_new`
    :param id: unique indentifer for motif
    :param r_struct: an RNAStructure instance to setup from

    :type structure: secondary_structure.Structure
    :type basepairs: list of secondary_structure.Basepairs
    :type ends: list of secondary_structure.Basepairs
    :type name: str
    :type path: str
    :type score: float
    :type end_ids: list of strs
    :type id: uuid.uuid1
    :type r_structure: secondary_structure.RNAStructure

    :attributes:

    `structure`: secondary_structure.Structure
        structure containing residue and chain information for
        this Motif instance
    `basepairs`: list of secondary_structure.Basepairs
        Basepairs between residues
    `ends`: list of secondary_structure.Basepairs
        Basepair ends where RNA structures can be connected
    `name`: str
        the name of the Pose
    `path`: str
        location of where Pose originated from, this is just
        a place holder for converting from pose.Pose
    `score` : float
        the score generated by motif_scorer.MotifScorer, estimates secondary
        structure stability
    `end_ids`: list of strs
        strings indenifying the secondary structure and sequence
        in the perspective of a given basepair end see
        :func:`assign_end_id_new`
    `motifs`: list of secondary_structure.Motifs
        all motifs that are contained in this Pose, excluding helices which
        need to be built and are stored in helices. Instead motifs stores all
        basepair steps as seperate motifs
    `helices`: list of secondary_structure.Motifs
        helical motifs that are more then 2 basepairs. needs to be build with
        :func:`Pose.build_helices`

    """

    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_end_ids",
        "_motifs"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, motifs):
        super(RNAStructure, self).__init__(structure, basepairs, ends, end_ids)
        self._motifs = motifs

    def __repr__(self):
        return "<secondary_structure.Pose( " + self.sequence() + " " + self.dot_bracket() + " )"

    def iter_motifs(self):
        return self._motifs.__iter__()

    def get_motif(self, m_uuid):
        """
        gets a motif by its unique indentifer

        :param m_id: motifs unique indentifier, Motif.id
        :type m_id: uuid.uuid1

        :returns: motif matching to its unique indentifier
        :rtype: secondary_structure.Motif
        """

        for m in self._motifs:
            if m.uuid == m_uuid:
                return m
        return None

    def replace_sequence(self, seq):
        """
        changes the sequence of structure.

        :param seq: the new sequence of structure
        :type seq: str

        :returns: None
        """

        super(self.__class__, self).replace_sequence(seq)

        for m in self._motifs:
            m.update()

    def update_motif(self, m_id):
        m = self.motif(m_id)
        for i, end in enumerate(m.ends):
            m.end_ids[i] = primitives.rna_structure.assign_end_id(m, end)

    def update(self):
        super(self.__class__, self).update()
        for m in self._motifs:
            m.update()

    def copy(self):
        """
        creates a deep copy of instance

        :returns: deep copy
        :rtype: secondary_structure.Pose
        """

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
        """
        generates a stringified verision of this instance.

        :returns: stringified verision of pose
        :rtype: str
        """

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
        """
        finds all basepair step motifs and builds helices from them. This is
        only required for external applications such as rosetta. helices get
        stored in member variable helices.
        """

        steps = []
        for m in self.motifs:
            if m.mtype != motif_type.HELIX:
                continue
            steps.append(m)

        seen = {}
        current = None
        helix_motifs = []
        found = 0
        while 1:
            helix_motifs = []
            current = None
            for m1 in steps:
                found = 0
                if m1 in seen:
                    continue
                for m2 in steps:
                    if m2 in seen:
                        continue
                    if m1 == m2:
                        continue
                    for end in m2.ends:
                        if m1.ends[0] == end:
                            found = 1
                            break
                if not found:
                    current = m1
                    break

            if found or current is None:
                break

            found = 1
            while found:
                seen[current] = 1
                helix_motifs.append(current)
                found = 0
                for m in steps:
                    if m in seen:
                        continue
                    if m.ends[0] == current.ends[1]:
                        current = m
                        found = 1
                        break

            res1, res2 = [], []
            bps, ends = [], []

            res1.append(helix_motifs[0].chains()[0].first())
            res2.append(helix_motifs[0].chains()[1].last())
            bps.append(helix_motifs[0].ends[0])
            ends.append(helix_motifs[0].ends[0])
            ends.append(helix_motifs[-1].ends[1])

            for m in helix_motifs:
                res1.append(m.chains()[0].last())
                res2.append(m.chains()[1].first())
                bps.append(m.ends[1])

            res2.reverse()
            chains = [Chain(res1), Chain(res2)]
            struc = Structure(chains)
            m = Motif(struc, bps, ends)
            self.helices.append(m)


class MotifDirectedGraph(motif_type_directed_graph.MotifTypeDirectedGraph):
    def __init__(self):
        super(self.__class__, self).__init__()
        self._sterics = 0
        self._motif_class_type = Motif
        self._motif_aligner = Aligner


def str_to_pose(s):
    """
    converts a pose from string generated from :func:`Pose.to_str`

    :param s: string created by Pose.to_str()
    :type s: str

    :return: pose from str
    :rtype: secondary_structure.Pose
    """

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

