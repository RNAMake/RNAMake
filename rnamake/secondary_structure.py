import uuid
import motif_type
import exceptions


class Residue(object):
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

    __slots__= ["name", "dot_bracket", "num", "chain_id", "uuid", "i_code"]

    def __init__(self, name, dot_bracket, num, chain_id, uuid, i_code=""):
        self.name, self.dot_bracket, self.num = name, dot_bracket, num
        self.uuid = uuid
        self.chain_id, self.i_code = chain_id, i_code

    def __repr__(self):
        return "<SecondaryStructureResidue('%s%d%s chain %s')>" % (
            self.name, self.num, self.i_code, self.chain_id)

    def to_str(self):
        """
        stringify residue object. can be converted back with
        :func:`str_to_residue`

        :return: stringifed verision of residue
        :rtype: str
        """

        return self.name + "," + self.dot_bracket + "," + str(self.num) + "," + \
               str(self.chain_id) + "," + str(self.i_code)

    def copy(self):
        """
        creates copy of current residue

        :return: copy of instatnce
        :rtype: secondary_structure.Residue
        """

        return Residue(self.name, self.dot_bracket, self.num, self.chain_id,
                       self.uuid, self.i_code)


class Chain(object):
    """
    secondary structure chain container. Contains a chain of connected
    residues. Chain should be from 5' to 3'.

    :param residues: the residues that are to be included in chain. Optional.
    :type residues: list of Residues.

    :attributes:
    `residues` : list of Residues
        Residues contained in current chain

    """

    __slots__ = ["residues"]

    def __init__(self, residues=None):
        self.residues = residues
        if self.residues is None:
            self.residues = []

    def __len__(self):
        return len(self.residues)

    def __repr__(self):
        seq = ""
        for r in self.residues:
            seq += r.name

        return "<SecondaryStructureChain( " + seq + ")"

    def __iter__(self):
        return self.residues.__iter__()

    def first(self):
        """
        gets the first residue in the chain

        :return: 5' end of chain
        :rtype: secondary_structure.Residue

        :examples:

        ..  code-block:: python

            >>> from rnamake.unittests import instances
            >>> c = instances.secondary_structure_chain()
            >>> c.first()
            <SecondaryStructureResidue('G13 chain A')>

        """

        if len(self.residues) == 0:
            raise exceptions.SecondaryStructureException(
                "cannot call first there are no residues in chain")

        return self.residues[0]

    def last(self):
        """
        gets the first residue in the chain

        :return: 5' end of chain
        :rtype: secondary_structure.Residue

        :examples:

        ..  code-block:: python

            >>> from rnamake.unittests import instances
            >>> c = instances.secondary_structure_chain()
            >>> c.last()
            <SecondaryStructureResidue('G24 chain A')>

        """

        if len(self.residues) == 0:
            raise exceptions.SecondaryStructureException(
                "cannot call last there are no residues in chain")

        return self.residues[-1]

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
        for r in self.residues:
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
        for r in self.residues:
            db += r.dot_bracket
        return db

    def to_str(self):
        """
        stringify chain object. can be converted back with
        :func:`str_to_chain`

        :return: stringifed verision of chain
        :rtype: str
        """

        s = ""
        for r in self.residues:
            s += r.to_str() + ";"
        return s

    def copy(self):
        """
        creates deep copy of chain instance

        :return: copy of chain
        :rtype: secondary_structure.Chain

        """

        residues = []
        for r in self.residues:
            residues.append(r.copy())
        return Chain(residues)


class Basepair(object):
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

    __slots__ = ["res1", "res2", "uuid"]

    def __init__(self, res1, res2, bp_uuid=None):
        self.res1, self.res2 = res1, res2
        self.uuid = bp_uuid
        if self.uuid is None:
            self.uuid = uuid.uuid1()

    def __repr__(self):
        return "<SecondaryStructureBasepair("+self.name()+")>"

    def name(self):
        """
        get name of basepair: which is the combined name of both residues
        seperated by a "-". The residue with the lower res number should
        come first

        :return: name of basepair
        :rtype: str

        :examples:

        ..  code-block:: python

            # build basepair from stratch
            >>> from rnamake.unittests import instances
            >>> b = instances.secondary_structure_basepair()
            >>> print b.res1
            <SecondaryStructureResidue('C10 chain A')>

            >>> print b.res2
            <SecondaryStructureResidue('G15 chain A')>

            >>> print b.name()
            A10-A15
        """

        str1 = self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code)
        str2 = self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

        if str1 < str2:
            return str1+"-"+str2
        else:
            return str2+"-"+str1

    def partner(self, r):
        """
        get the other basepairing partner of a residue will throw an error
        if the supplied residue is not contained within this basepair

        :param res: the residue that you want to get the partner of
        :type res: secondary_structure.Residue object
        """

        if   r == self.res1:
            return self.res2
        elif r == self.res2:
            return self.res1
        else:
            raise exceptions.SecondaryStructureException(
                "call partner with a residue not in basepair")


class Structure(object):
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

    __slots__ = ["chains"]

    def __init__(self, chains=None, sequence="", dot_bracket=""):
        self.chains = []
        if chains is not None:
            self.chains = chains

        if len(sequence) != 0 and len(dot_bracket) != 0:
           self.chains = self._setup_chains(sequence, dot_bracket)

    def _setup_chains(self, sequence, dot_bracket):
        """
        setup function for turning a string sequence and secondary structure
        into a structure object.

        :param sequence: sequence of RNA
        :param dot_bracket: dot bracket of RNA

        :type sequence: str
        :type dot_bracket: str
        """

        chains = []
        residues = []

        if len(dot_bracket) != len(sequence):
            raise exceptions.SecondaryStructureException(
                "sequence and dot bracket are not the same length")

        if dot_bracket[0] != '(' and dot_bracket[0] != '.' and dot_bracket != '&':
            raise exceptions.SecondaryStructureException(
                "secondary structure is not valid did you flip seq and ss?")

        count = 1
        chains_ids = "ABCDEFGHIJKLMNOPQRSTUVWXZ"
        valid_seq = "AGUCTN&+-"
        chain_i = 0
        for i in range(len(sequence)):
            if sequence[i] not in valid_seq:
                raise exceptions.SecondaryStructureException(
                    sequence[i] + " is not a valid secondary_structure element")

            if sequence[i] != "&" and sequence[i] != "+" and sequence[i] != "-":
                r = Residue(sequence[i], dot_bracket[i], count,
                            chains_ids[chain_i], uuid.uuid1())
                residues.append(r)
                count += 1
            else:
                chain_i += 1
                chains.append(Chain(residues))
                # unlikely but hit max chains
                if chain_i == len(chains_ids)-1:
                    chain_i = 0
                residues = []

        if len(residues) > 0:
            chains.append(Chain(residues))

        return chains

    def residues(self):
        """
        Concats all residue objects from all Chain objects intos a unified
        list to be able to easily iterate through.

        :return: List of secondary_structure.Residue objects
        """
        res = []
        for c in self.chains:
            res.extend(c.residues)
        return res

    def sequence(self):
        """
        Concats the sequence of each Chain into one sequence for the entire
        RNA

        :return: sequence of structure
        :rtype: seq
        """
        sequences = [x.sequence() for x in self.chains]
        return "&".join(sequences)

    def dot_bracket(self):
        """
        Concats the secondary structure in the form of dot bracket notation
        of each Chain into one sequence for the entire RNA

        :return: sequence of structure
        :rtype: seq
        """
        dot_brackets = [x.dot_bracket() for x in self.chains]
        return "&".join(dot_brackets)

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        """
        find a residue based on residue num, chain_id, insert_code and uuid
        will return an error if more then one residue matches search to avoid
        confusion. Will return None is nothing matches search.

        :param num: residue number
        :param chain_id: what chain the residue belongs to
        :param i_code: the insertation code of the residue
        :param uuid: the unique indentifier that each residue is given

        :type num: int
        :type chain_id: str
        :type i_code: str
        :type uuid: uuid

        :return: Residue object
        :rtype: residue.Residue

        :examples:

        .. code-block:: python

            >>> from rnamake import secondary_structure
            >>> s = secondary_structure.Structure(sequence="GCGAAAACGC",
                                                  dot_bracket="(((....)))")
            >>> s.get_residue(num=1)
            <SecondaryStructureResidue('G1 chain A')>
        """

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
        """
        creates a deep copy of structure instance

        :return: copy of structure
        :rtype: secondary_structure.Structure
        """

        new_chains = [ c.copy() for c in self.chains]
        return Structure(chains=new_chains)

    def to_str(self):
        """
        generates a stringified verision of this instance.

        :returns: stringified verision of structure
        :rtype: str
        """
        s = ""
        for c in self.chains:
            s += c.to_str() + "|"
        return s


class RNAStructure(object):
    """
    Complete secondary structure container for representing a RNA. Contains
    both the sequence indentity of each residue with its corresponding dot
    bracket notation symbol but also includes basepair objects to represent the
    pairs between residues. This class parallels the
    rna_structure.RNAStructure class for describing RNA with 3D coordinates.
    Having a parallel class for secondary structure makes it simple to move
    between secondary structure and full atom representations of RNA.

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



    """

    def __init__(self, structure=None, basepairs=None, ends=None,
                 name="assembled", path="assembled", score=0, end_ids=None):
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
        self.helices = []

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


def str_to_residue(s):
    """
    converts a residue from string generated from Resiude.to_str()

    :param s: string created by Residue.to_str()
    :type s: str

    :return: residue from str
    :rtype: secondary_structure.Residue
    """

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

