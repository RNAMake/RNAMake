import numpy as np
import uuid

# RNAMake imports
from rnamake import transform
from rnamake import util
from rnamake import settings
from rnamake import basic_io
from rnamake import secondary_structure
from rnamake import bead
from rnamake import motif_state
from rnamake import primitives
from rnamake import exceptions
from rnamake import x3dna
from rnamake import user_warnings


class Atom(object):
    """Stores atomic information from pdb file, design is to be extremely
    lightweight only storing the atom name and coordinates.

    :param name: name of atom
    :param coords: 3d coordinates of atom's position

    :type name: str
    :type coords: numpy.array


    :attributes:

    `__name` : str
        Atomic name
    `__coords` : np.array
        Atomic coordinates

    :examples:

    .. code-block:: python

        >>> a = Atom("P",[1.0,2.0,3.0])
        >>> a.name
        P

        >>> a.coords
        [1.0 2.0 3.0]

        >>> print a
        <Atom(name ='P', coords='1.0 2.0 3.0')>

    """
    __slots__ = ["__name", "__coords"]

    def __init__(self, name, coords):
        """returns new atom.Atom object"""

        # stop from users giving a list instead, this will mess stuff up in
        # unexpected ways
        if type(coords) is not np.ndarray:
            raise exceptions.AtomException(
                "invalid creation of atom, coords needs to be np.array")

        self.__name, self.__coords = name, coords

    @classmethod
    def from_str(cls, s):
        """
        converts string to atom.Atom object format "AtomName X Y Z"

        :params s: string containing atom elements
        :type s: str

        .. code-block:: python

            >>> Atom.from_str("P 1.0 2.0 3.0")
            <Atom(name='P', coords='1.0 2.0 3.0')>
        """

        spl = s.split()
        if len(spl) != 4:
            raise exceptions.AtomException(
                "did not get correct number of elements in string to build atom")
        coords = [float(x) for x in spl[1:]]
        return cls(spl[0], np.array(coords))

    @classmethod
    def copy(cls, a):
        """Deep copies the an atom instance.

        :returns: an Atom object

        :examples:

        .. code-block:: python

            >>> a = Atom("P",[1.0,2.0,3.0])
            >>> a_copy = Atom.copy(a)
            >>> print a_copy.name
            P

        """

        return cls(a.__name, np.array(a.__coords, copy=True))

    def __repr__(self):
        """returns string representation of object"""

        coords = basic_io.point_to_str(self.coords)
        return "<Atom(name='%s', coords='%s')>" % (self.name, coords)

    def transform(self, t):
        self.__coords = np.dot(self.__coords, t.rotation().T) + t.translation()

    def move(self, p):
        self.__coords += p

    # getters
    def get_copy(self):
        return Atom.copy(self)

    def get_str(self):
        """returns string version of atom.

        :returns: str

        :examples:

        .. code-block:: python

            >>> atom = atom.Atom("H1", np.array([0, 1, 2]))
            >>> string = atom.to_str()
            "H1 0.0 1.0 2.0"
        """
        return self.__name + " " + basic_io.point_to_str(self.__coords)

    def get_pdb_str(self, acount=1):
        """
        prints the current atom to a pdb string

        :param acount: the atom number of the atom in pdb file
        :type acount: int

        :examples:

        .. code-block:: python

            >>> a = Atom("P",[1.0,2.0,3.0])
            >>> a.to_pdb_str()
            ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P

            >>> a.to_pdb_str(10)
            ATOM     10  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P
        """
        if self is None:
            return ""

        s = "ATOM {:6d}  P   C   A   1 {:11.3f}{:8.3f}{:8.3f}  1.00 62.18           P\n".format(
            acount,
            self.__coords[0],
            self.__coords[1],
            self.__coords[2])
        return s

    def get_coords(self):
        return np.copy(self.__coords)

    def get_name(self):
        return self.__name


class Residue(primitives.Residue):
    """
    Store residue information from pdb file, stores all Atom objects that
    belong to residue. Implementation is designed to be extremely lightweight.

    :param rtype: residue type, stores information about residue
    :type rtype: residue_type.ResidueType

    :param name: residue name
    :type name: str

    :param num: residue num
    :type num: int

    :param chain_id: chain identification
    :type chain_id: str

    :param i_code: insertion code, optional argument if not supplied will be
        set to
    :type i_code: str

    :attributes:

    `_atoms` : atom.Atom
        holds all atoms that belong to this residue object
    `_name` : str
        name of residue, ex. ADE, GUA etc
    `_num` : int
        residue num
    `_rtype` : residue_type.ResidueType
        Information about residue type each nucleic acid has its own type
    `_chain_id` : str
        chain indentification string, ex. 'A' or 'B'
    `_i_code`: str
        residue insertion code
    `_uuid` : uuid.uuid1
        the unique id to indentify residue

    :examples:

    .. code-block:: python

        # generating a new residue
        >>> from all_atom import Residue
        >>> rts = residue_type.ResidueTypeSet()
        >>> rtype = rts.get_rtype_by_resname("ADE")
        >>> r = Residue(rtype, "ADE", 1, "A")
        >>> print r.name
        A

        #using test residue
        >>> import rnamake.unittests.instances
        >>> r = rnamake.unittests.instances.residue()
        >>> print r.name
        G

        >>> a = r.get_atom("C1'")
        >>> print a.coords
        [-23.806 -50.289  86.732]

        >>> r.get_beads()
        [<Bead(btype='SUGAR', center='-24.027 -48.5001111111 86.368')>, <Bead(btype='BASE', center='-21.2186363636 -52.048 85.1157272727')>]

        #a fast way of saving coordinate information to file
        >>> r.to_str()
        "GUA,G,103,A,,N,N,N,O5' -26.469 -47.756 84.669,C5' -25.05 -47.579 84.775,C4' -24.521 -48.156 86.068,O4' -24.861 -49.568 86.118,C3' -23.009 -48.119 86.281,O3' -22.548 -46.872 86.808,C1' -23.806 -50.289 86.732,C2' -22.812 -49.259 87.269,O2' -23.167 -48.903 88.592,N1 -19.538 -52.485 85.025,C2 -19.717 -51.643 86.097,N2 -18.624 -51.354 86.809,N3 -20.884 -51.124 86.445,C4 -21.881 -51.521 85.623,C5 -21.811 -52.356 84.527,C6 -20.546 -52.91 84.164,O6 -20.273 -53.677 83.228,N7 -23.063 -52.513 83.947,C8 -23.858 -51.786 84.686,N9 -23.21 -51.159 85.722,"

        #get PDB formmated coordinates back out
        >>> r.to_pdb_str()
        ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00
        ATOM      2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00
        ATOM      3 C4'  G   A 103     -24.521 -48.156  86.068  1.00  0.00
        ATOM      4 O4'  G   A 103     -24.861 -49.568  86.118  1.00  0.00
        ATOM      5 C3'  G   A 103     -23.009 -48.119  86.281  1.00  0.00
        ATOM      6 O3'  G   A 103     -22.548 -46.872  86.808  1.00  0.00
        .
        .
        .
    """

    __slots__ = [
        "_atoms",
        "_rtype",
        "_name",
        "_num",
        "_chain_id",
        "_i_code",
        "_uuid",
        "_beads"]

    def __init__(self, atoms, rtype, name, num, chain_id, i_code=None, r_uuid=None):
        self._rtype = rtype
        self._atoms = []
        self.__setup_atoms(atoms)
        self._beads = []

        super(self.__class__, self).__init__(name, num, chain_id,
                                             i_code=i_code, r_uuid=r_uuid)

    @classmethod
    def from_str(cls, s, rts):
        """
        creates an residue from string generated from
        :func:`rnamake.residue.Residue.to_str`

        :param s: string containing stringifed residue
        :type s: str

        :returns: unstringifed residue object
        :rtype: residue.Residue

        """

        spl = s.split(",")
        atoms = []
        for i in range(5, len(spl)-1):
            if spl[i] == "N":
                a = None
            else:
                a = Atom.from_str(spl[i])
            atoms.append(a)

        rtype = rts.get_type(spl[0])
        r = cls(atoms, rtype, spl[1], int(spl[2]), spl[3], spl[4])
        return r

    @classmethod
    def copy(cls, r, new_uuid=0, build_beads=1, given_uuid=None):
        """
        performs a deep copy of Residue object

        :rtype: Residue

        :examples:

        .. code-block:: python

            >>> import rnamake.unittests.instances
            >>> r = rnamake.unittests.instances.residue()
            >>> r_copy = Residue.copy(r)
            >>> r_copy.name
            G

        """

        atoms = []
        for a in r:
            if a is not None:
                atoms.append(Atom.copy(a))

        r_uuid = r._uuid
        if new_uuid:
            r_uuid = uuid.uuid1()
        if given_uuid:
            r_uuid = given_uuid

        new_r = cls(atoms, r._rtype, r._name, r._num, r._chain_id,
                    r._i_code, r_uuid)

        if build_beads:
            if r.get_num_beads() > 0:
                new_r.build_beads()

        return new_r

    def __repr__(self):
        return "<Residue('%s%d%s chain %s')>" % (
            self._name, self._num, self._i_code, self._chain_id)

    # iterators
    def __iter__(self):
        return self._atoms.__iter__()

    def iter_beads(self):
        return self._beads.__iter__()

    # modifiers
    def build_beads(self):
        self._beads = self.__get_beads()

    def has_atom(self, atom_name=None, index=None):
        pos = None
        if atom_name is not None:
            pos = self.__get_atom_position(atom_name)
            if pos is None:
                return False

        if index is not None:
            if index >= len(self._atoms):
                return False
            pos = index

        if self._atoms[pos] is None:
            return False
        else:
            return True

    def move(self, p):
        for a in self._atoms:
            if a is not None:
                a.move(p)
        for b in self._beads:
            b.move(p)

    def transform(self, t):
        for a in self._atoms:
            if a is not None:
                a.transform(t)
        for b in self._beads:
            b.transform(t)

    def remove_beads(self):
        self._beads = []

    # getters
    def get_atom(self, atom_name=None, index=None):
        """
        get atom object by its name

        :param atom_name: name of atom
        :type atom_name: str

        Examples:

        .. code-block:: python

            >>> r = rnamake.unittests.instances.residue()
            >>> a = r.get_atom("C1'")
            >>> print a.coords
            [-23.806 -50.289  86.732]
        """
        if atom_name is not None:
            pos = self.__get_atom_position(atom_name)
            if pos is None:
                raise exceptions.ResidueException(
                    "atom: " + atom_name + " does not exist in residue")

        if index is not None:
            if index >= len(self._atoms):
                raise exceptions.ResidueException(
                    "atom pos: " + str(index) + " does not exist in residue")
            pos = index

        if self._atoms[pos] is None:
            if atom_name is not None:
                raise exceptions.ResidueException(
                    "atom: " + atom_name + " is not initialized in residue")
            else:
                raise exceptions.ResidueException(
                    "atom: " + str(pos) + " is not initialized in residue")

        return self._atoms[pos]

    def get_bead(self, btype):
        for b in self._beads:
            if b.btype == btype:
                return b
        return None

    def get_center(self):
        return util.center(self._atoms)

    def get_coords(self, atom_name=None, index=None):
        return self.get_atom(atom_name, index).get_coords()

    def get_copy(self, new_uuid=0, build_beads=1, given_uuid=None):
        return Residue.copy(self, new_uuid, build_beads, given_uuid)

    def get_pdb_str(self, acount=1, return_acount=0, rnum=-1, chain_id=""):
        """
        returns pdb formatted string of residue's coordinate information

        :param acount: current atom index, default: 1
        :param return_acount: final atom index after current atoms, default: 0
        :param rnum: starting residue number, default: -1
        :param chain_id: the chain id of the chain, i.e. "A", "B" etc

        :type  acount: int
        :type  return_acount: int
        :type  rnum: int
        :type  chain_id: str

        :rtype: str

        :examples:

        .. code-block:: python

            >>> import rnamake.unittests.instances
            >>> r = rnamake.unittests.instances.residue()
            >>> r.to_pdb_str()
            ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00
            ATOM      2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00
            ATOM      3 C4'  G   A 103     -24.521 -48.156  86.068  1.00  0.00
            ATOM      4 O4'  G   A 103     -24.861 -49.568  86.118  1.00  0.00
            ATOM      5 C3'  G   A 103     -23.009 -48.119  86.281  1.00  0.00
            ATOM      6 O3'  G   A 103     -22.548 -46.872  86.808  1.00  0.00
            .
            .
            .
        """

        num = self._num
        cid = self._chain_id
        if rnum != -1:
            num = rnum
        if chain_id != "":
            cid = chain_id

        s = ""
        for a in self._atoms:
            if a is None:
                continue
            s += basic_io.PDBLINE_GE100K % \
                 ('ATOM', acount, a.name, '', self.get_short_name(), cid,
                  num, '', a.coords[0], a.coords[1], a.coords[2], 1.00,
                  0.00, '', '')
            acount += 1

        if return_acount:
            return s, acount
        else:
            return s

    def get_pdb(self, fname="residue.pdb"):
        """
        Writes a PDB string formmated verision of this Residue object to file

        :param fname: filename of output PDB file, default="residue.pdb"

        :type  fname: str
        :return: None
        """

        f = open(fname, "w")
        s = self.to_pdb_str()
        f.write(s)
        f.close()

    def get_num_beads(self):
        return len(self._beads)

    def get_num_atoms(self):
        count = 0
        for a in self._atoms:
            if a is not None:
                count += 1
        return count

    def get_short_name(self):
        """gets letter of residue, i.e. A or G etc

        :return: letter for residue
        :rtype: str
        """
        return self._rtype.short_name

    def get_state(self):
        beads = [bead.Bead.copy(b) for b in self._beads]
        return motif_state.Residue(self._name, self._num, self._chain_id,
                                   self._i_code, beads, self._uuid)

    def get_str(self):
        """
        stringifes residue object

        :returns: stringified residue object

        .. code-block:: python

            >>> import rnamake.unittests.instances
            >>> r = rnamake.unittests.instances.residue()
            >>> r.to_str()
            "GUA,G,103,A,,N,N,N,O5' -26.469 -47.756 84.669,C5' -25.05 -47.579 84.775,C4' -24.521 -48.156 86.068,O4' -24.861 -49.568 86.118,C3' -23.009 -48.119 86.281,O3' -22.548 -46.872 86.808,C1' -23.806 -50.289 86.732,C2' -22.812 -49.259 87.269,O2' -23.167 -48.903 88.592,N1 -19.538 -52.485 85.025,C2 -19.717 -51.643 86.097,N2 -18.624 -51.354 86.809,N3 -20.884 -51.124 86.445,C4 -21.881 -51.521 85.623,C5 -21.811 -52.356 84.527,C6 -20.546 -52.91 84.164,O6 -20.273 -53.677 83.228,N7 -23.063 -52.513 83.947,C8 -23.858 -51.786 84.686,N9 -23.21 -51.159 85.722,"
        """
        s =  self._rtype.name + "," + self._name + "," + str(self._num) + ","
        s += self._chain_id + "," + self._i_code + ","
        for a in self._atoms:
            if a is None:
                s += "N,"
            else:
                s += a.get_str() + ","
        return s

    def get_sugar_bead(self):
        for b in self._beads:
            if b.btype == bead.BeadType.SUGAR:
                return b
        return None

    # private
    def __get_atom_position(self, atom_name):
        if self._rtype.is_valid_atom(atom_name):
            return self._rtype.atom_index(atom_name)
        else:
            return None

    def __setup_atoms(self, atoms):
        """
        put atoms in correct positon in internal atom list, also corrects some
        named atom names to their correct name

        :param atoms: list of atom objects that are to be part of this residue
        :type atoms: list of Atom objects
        """

        self._atoms = [None for x in range(0, len(self._rtype))]
        for a in atoms:
            if a is None:
                continue
            name_change = self._rtype.get_correct_atom_name(a)
            if name_change is not None:
                a = Atom(name_change[1], a.get_coords())
            if self._rtype.is_valid_atom(a.get_name()):
                pos = self._rtype.atom_index(a.get_name())
                self._atoms[pos] = a

    def __get_beads(self):
        """
        Generates steric beads required for checking for steric clashes between
        motifs. Each residues has three beads modeled after the typical three
        bead models used in coarse grain modeling. The three beads are:

        Phosphate:  P, OP1, OP2\n
        Sugar    :  O5',C5',C4',O4',C3',O3',C1',C2',O2'\n
        Base     :  All remaining atoms

        if there are for example no phosphate atoms only 2 beads will be returned.

        .. code-block:: python

            >>> import rnamake.unittests.instances
            >>> r = rnamake.unittests.instances.residue()
            >>> r.get_beads()
            [<Bead(btype='SUGAR', center='-24.027 -48.5001111111 86.368')>, <Bead(btype='BASE', center='-21.2186363636 -52.048 85.1157272727')>]

        """
        phos_atoms, sugar_atoms, base_atoms = [], [], []

        for i, a in enumerate(self._atoms):
            if a is None:
                continue
            if i < 3:
                phos_atoms.append(a)
            elif i < 12:
                sugar_atoms.append(a)
            else:
                base_atoms.append(a)

        beads = []
        types = [bead.BeadType.PHOS, bead.BeadType.SUGAR, bead.BeadType.BASE]
        for i, alist in enumerate([phos_atoms, sugar_atoms, base_atoms]):
            if len(alist) > 0:
                beads.append(bead.Bead(util.center(alist), types[i]))

        return beads


class Chain(primitives.Chain):
    """
    Stored chain information from pdb file. Stores all residues in chain.
    Implementation is designed to be extremely lightweight. To connect residues
    into chains it highly adviced that you use
    :func:`connect_residues_into_chains`

    :param residues: the residues that are to be included in this chain
    :type residues: list of residue.Residue objects

    :attributes:

    `residues` : List of Residue objects
        The list of residues that belong to this chain will always be in
        5' to 3' order

    :examples:

    ..  code-block:: python

        >>> import rnamake.unittests.instances
        >>> c = rnamake.unittests.instances.chain()
        >>> c.first()
        <Residue('G103 chain A')>

        >>> c.last()
        <Residue('C260 chain A')>

        >>> len(c)
        157

        >>> cs = c.subchain(1, 10)
        >>> len(cs)
        9

        >>> cs.first()
        <Residue('A104 chain A')>

        >>> cs2 = c.subchain(start_res=c.residues[10], end_res=c.residues[15])
        >>> len(cs2)
        6

        >>> c.to_pdb_str()
        ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00
        ATOM      2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00
        ATOM      3 C4'  G   A 103     -24.521 -48.156  86.068  1.00  0.00
        ATOM      4 O4'  G   A 103     -24.861 -49.568  86.118  1.00  0.00
        ATOM      5 C3'  G   A 103     -23.009 -48.119  86.281  1.00  0.00
        ATOM      6 O3'  G   A 103     -22.548 -46.872  86.808  1.00  0.00
        ATOM      7 C1'  G   A 103     -23.806 -50.289  86.732  1.00  0.00
        ATOM      8 C2'  G   A 103     -22.812 -49.259  87.269  1.00  0.00
        .
        .
        .

    """
    __slots__ = ["_residues"]

    def __repr__(self):
        if len(self._residues) == 0:
            return "<Chain( First: None\n\t  Last:  None\n\t  Size: 0)>"

        return "<Chain( First: %s\n\t  Last:  %s\n\t  Size: %s)>" %\
            (self.first(), self.last(), len(self._residues))

    @classmethod
    def from_str(cls, s, rts):
        """
        creates an chain from string generated from
        :func:`rnamake.chain.Chain.to_str`

        :param s: string containing stringifed chain
        :type s: str

        :returns: unstringifed chain object
        :rtype: chain.Chain
        """
        spl = s.split(";")
        residues = []
        for r_str in spl[:-1]:
            r = Residue.from_str(r_str, rts)
            residues.append(r)
        return cls(residues)

    @classmethod
    def copy(cls, c, new_uuid=0):
        """
        Creates a deepcopy of the this chain object.

        :return: Chain object
        """
        residues = [Residue.copy(r, new_uuid) for r in c]
        return cls(residues)

    def move(self, p):
        for r in self._residues:
            r.move(p)

    def transform(self, t):
        for r in self._residues:
            r.transform(t)

    # getters
    def get_copy(self, new_uuid=0):
        return Chain.copy(self, new_uuid)

    def get_pdb_str(self, acount=1, return_acount=0, rnum=-1, chain_id=""):
        """
        creates a PDB string formatted verision of this Chain object.

        :param acount: current atom index, default: 1
        :param return_acount: final atom index after current atoms, default: 0
        :param rnum: starting residue number, default: -1
        :param chain_id: the chain id of the chain, i.e. "A", "B" etc

        :type  acount: int
        :type  return_acount: int
        :type  rnum: int
        :type  chain_id: str

        :return: str
        """

        s = ""
        for r in self._residues:
            r_str, acount = r.get_pdb_str(acount, 1, rnum, chain_id)
            if rnum != -1:
                rnum += 1
            s += r_str

        # TODO fix returns should only be one possibility
        if return_acount:
            return s, acount
        else:
            return s

    def get_pdb(self, fname="chain.pdb", rnum=-1, chain_id=""):
        """
        Writes a PDB string formmated verision of this Chain object to file

        :param fname: filename of output PDB file, default="chain.pdb"
        :param rnum: starting residue number, default: -1
        :param chain_id: the chain id of the chain, i.e. "A", "B" etc

        :type  fname: str
        :type  rnum: int
        :type  chain_id: str

        :return: None
        """
        f = open(fname, "w")
        f.write(self.get_pdb_str(rnum=rnum, chain_id=chain_id))
        f.close()

    def get_state(self):
        residues = [r.get_state() for r in self._residues]
        return motif_state.Chain(residues)

    def get_str(self):
        """
        Stringifes Chain object

        :return: str
        """

        s = ""
        for r in self._residues:
            s += r.get_str() + ";"
        return s


class Structure(primitives.Structure):
    """Stores 3D structure information from a pdb file. Stores all chains,
    residues and atoms objects. Implementation is designed to be extremely
    lightweight and capable of performing fast transformations. to load a PDB
    formated file into a Structure object use
    :func:`structure_from_pdb`

    :attributes:

    `residues_` : list of residue.Residues
        These residues belong to the current structure

    :examples:

    ..  code-block:: python

        # load structure from pdb formatted file
        >>> import rnamake.unittests.files
        >>> s = structure_from_pdb(rnamake.unittests.files.P4P6_PDB_PATH)

        # load structure from test instance
        >>> import rnamake.unittests.instances
        >>> s = rnamake.unittests.instances.structure()

        # get specific residue
        >>> r = s.get_residue(num=106)
        >>> print r
        <Residue('U106 chain A')>

        # get residue using its unique indentifer
        >>>s.get_residue(uuid=r.uuid)
        <Residue('U106 chain A')>

        >>> s.chains()
        [<Chain( First: <Residue('G103 chain A')>
                Last:  <Residue('C260 chain A')>
                Size: 157)>]

        >>> len(s.get_beads())
        470

        # exclude beads from first residue. This can be useful if you only
        # need sterics from a part of the structure
        >>> len(s.get_beads(excluded_res=[s.residues()[0]]))
        468

    """

    __slots__ = [
        "_residues",
        "_chain_cuts"]

    def __init__(self, residues, chain_cuts):
       super(self.__class__, self).__init__(residues, chain_cuts)

    @classmethod
    def from_str(cls, s, rts):
        """
        creates an structure from string generated from
        :func:`rnamake.structure.Structure.to_str`
        :param s: string containing stringifed structure
        :type s: str
        :returns: unstringifed structure object
        :rtype: structure.Structure
        """

        spl = s.split(";")
        residues = []
        for r_str in spl[:-2]:
            r = Residue.from_str(r_str, rts)
            residues.append(r)
        chain_cuts = [int(x) for x in spl[-2].split()]
        return cls(residues, chain_cuts)

    @classmethod
    def copy(cls, s, new_uuid=0):
        """
        creates a deep copy of this structure

        :returns: copy of Structure object
        :rtype: Structure
        """
        residues = []
        for r in s._residues:
            cr = Residue.copy(r, new_uuid)
            residues.append(cr)

        return cls(residues, s._chain_cuts)

    def __repr__(self):
        return """<Structure(#chains: %s, #residues: %s)>""" %\
               (len(self.get_chains()), len(self._residues))

    def move(self, p):
        """
        Moves atomic coordinates based on a specified translation. If you do
        not supply and np.array and instead give it a list it will not work!

        :param p: Amount to 3D space to displace each atom
        :type p: numpy.array

        :return: None

        :examples:

        .. code-block:: python

            >>> import rnamake.unittests.instances
            >>> import numpy as np

            # print out original and moved structure of P4-P6 domain
            >>> s = rnamake.unittests.instances.structure()
            >>> s.to_pdb("org.pdb")
            >>> s.move(np.array([50, 0, 0]))
            >>> s.to_pdb("moved.pdb")

        .. figure:: ../_static/img/move_test.png
            :align:   center

            Result of example above, Red original, Blue moved.

        """

        for r in self._residues:
            r.move(p)

    def transform(self, t):
        """
        Transforms atomic coordinates based on a rotation and translation

        :param t: Transformation to be performed on atoms
        :type t: transform.Transformation

        :return: None

        :examples:

        .. code-block:: python

            # load structure from test instance
            >>> import rnamake.unittests.instances
            >>> s = rnamake.unittests.instances.structure()

            >>> import rnamake.util
            >>> rnamake.util.center(s.atoms())
            array([ -8.49624784, -59.88166488,  80.1866473 ])

            # load indenity transform (no changes)
            >>> t = rnamake.unittests.instances.transform_indentity()
            >>> s.transform(t)
            >>> rnamake.util.center(s.atoms())
            array([ -8.49624784, -59.88166488,  80.1866473 ])

            # load random transform, write pdbs to disk
            >>> t = rnamake.unittests.instances.transform_random()
            >>> s.to_pdb("org.pdb")
            >>> s.transform(t)
            >>> s.to_pdb("transformed.pdb")

        .. figure:: ../_static/img/transform_test.png
            :align:   center

            Result of example above, Red original, Blue transformed.


        """
        for r in self._residues:
            r.transform(t)

    # getters
    def get_copy(self, new_uuid=0):
        return Structure.copy(self, new_uuid)

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

    def get_pdb_str(self, renumber=-1):
        """
        creates a PDB string formatted verision of this Structure object.

        :param renumber: what should the first residue be numbered. -1 is
            to NOT renumber, Default=-1.
        :type renumber: int

        :return: str
        """
        acount = 1
        s = ""
        c_names = "ABCDEFGHIJKLM"
        rnum = -1
        chain_id = ""
        if renumber != -1:
            chain_id = "A"
            rnum = 1

        for i, c in enumerate(self.get_chains()):
            c_str, acount = c.get_pdb_str(acount, 1, rnum, chain_id)
            if renumber != -1:
                #    chain_id = c_names[i+1]
                rnum += len(c)
            s += c_str
            s += "TER\n"
        return s

    def get_pdb(self, fname="structure.pdb", renumber=-1):
        """
        write structure to pdb file

        :param fname: name of the file of the pdb file you want to write to
        :param renumber: what should the first residue be numbered. -1 is
            to NOT renumber, Default=-1.

        :type fname: str
        :type renumber: int

        :return: None

        """
        f = open(fname, "w")
        f.write(self.to_pdb_str(renumber))
        f.close()

    def get_state(self):
        residues = [ r.get_state() for r in self._residues ]
        return motif_state.Structure(residues, self._chain_cuts)

    def get_str(self):
        """
        Stringifes Structure object

        :return: str
        """

        s = ""
        for r in self._residues:
            s += r.get_str() + ";"
        s += " ".join([str(x) for x in self._chain_cuts]) + ";"
        return s


class Basepair(primitives.Basepair):
    """
    :param res1: First residue in basepair
    :param res2: Second residue in basepair
    :param r: Reference frame of basepair
    :param bp_type: X3dna basepair type, default "c..."

    :type res1: residue.Residue
    :type res2: residue.Residue
    :type r: np.array
    :type bp_type: str

    :attributes:

    `res1` : residue.Residue
        First residue in basepair
    `res2` : residue.Residue
        Second residue in basepair
    `bp_type`: str
        X3dna basepair type
    `atoms`: list of atom.Atoms
        All atoms from both res1 and res2. Grouped together for easy manipulation
    `bp_state`: BasepairState
        Contains the orientation, origin and sugar. This is all information
        required to align to this basepair and can be used independently from
        the reset of the basepair.
    `uuid`: uuid.uuid1()
        unique id to indentify this basepair when locating it in a motif or
        pose

    :examples:

    ..  code-block:: python

        # build basepair from stratch
        >>> from rnamake.unittests import instances
        >>> from rnamake import basepair
        >>> import numpy as np
        >>> s = instances.structure()
        >>> b = basepair.Basepair(s.get_residue(num=103),
                                  s.get_residue(num=104),
                                  np.eye(3))
        >>> print b
        <Basepair(A103-A104)>

        # loading test basepair
        >>> b = instances.basepair()
        >>> print b
        <Basepair(A13-A12)>

        # primary axis of orientation, used to align to this basepair
        >>> print b.r()
        [[  1.00000001e+00   1.00000001e-04  -9.99999990e-09]
         [ -1.00000002e-04   1.00000001e+00  -9.99999990e-05]
         [ -1.99999995e-08   9.99999955e-05   1.00000000e+00]]

        # center of mass of basepair
        >>> print b.d()
        [ 0.1956032   0.69256601  0.0930465 ]

    """

    __slots__ = [
        "_res1_uuid",
        "_res2_uuid",
        "_r",
        "_d",
        "_sugars",
        "_name",
        "_bp_type",
        "_x3dna_bp_type",
        "_uuid"]

    def __init__(self, res1_uuid, res2_uuid, r, d, sugars, name,
                 x3dna_bp_type=None, bp_type=None, bp_uuid=None):
        self._res1_uuid, self._res2_uuid = res1_uuid, res2_uuid
        self._r = r
        self._d = d
        self._sugars = sugars
        self._name = name
        self._bp_type = bp_type
        self._uuid = bp_uuid
        self._x3dna_bp_type = x3dna_bp_type

        if self._bp_type is None:
            self._bp_type = primitives.basepair.BasepairType.NC

        if self._x3dna_bp_type is None:
            self._x3dna_bp_type = x3dna.X3dnaBPType.cDDD

        if self._uuid is None:
            self._uuid = uuid.uuid1()

    @classmethod
    def from_str(cls, s, res1_uuid, res2_uuid):
        spl = s.split(";")
        d = basic_io.str_to_point(spl[0])
        r = basic_io.str_to_matrix(spl[1])
        sugars = basic_io.str_to_points(spl[2])

        return cls(res1_uuid, res2_uuid, r, d, sugars, spl[3], int(spl[4]), int(spl[5]))

    @classmethod
    def copy(cls, bp):
        sugars = [np.copy(bp._sugars[0]), np.copy(bp._sugars[1])]
        return cls(bp._res1_uuid, bp._res2_uuid, np.copy(bp._r), np.copy(bp._d),
                   sugars, bp._name, bp._x3dna_bp_type, bp._bp_type, bp._uuid)

    @classmethod
    def copy_with_new_uuids(cls, bp, res1_uuid, res2_uuid, bp_uuid=None):
        if bp_uuid is None:
            bp_uuid = uuid.uuid1()
        sugars = [np.copy(bp._sugars[0]), np.copy(bp._sugars[1])]
        return cls(res1_uuid, res2_uuid, np.copy(bp._r), np.copy(bp._d),
                   sugars, bp._name, bp._x3dna_bp_type, bp._bp_type, bp_uuid=bp_uuid)

    def __repr__(self):
          return "<Basepair("+self._name + ")>"

    def diff(self, bp):
        diff = util.distance(self.d, bp.d)
        diff += self._rot_diff(bp) * 2
        return diff

    # modifiers
    def flip_res(self):
        self._res1_uuid, self._res2_uuid = self._res2_uuid, self._res1_uuid

    def flip(self):
        self.r[1] = -self.r[1]
        self.r[2] = -self.r[2]

    def transform(self, t):
        r_T = t.rotation().T

        new_r = util.unitarize(np.dot(self._r, r_T))
        new_sugars = []
        for s in self._sugars:
            new_sugars.append(np.dot(s, r_T) + t.translation())
        new_origin = np.dot(self._d, r_T) + t.translation()

        self._r = new_r
        self._d = new_origin
        self._sugars = new_sugars

    def move(self, p):
        self._sugars[0] += p
        self._sugars[1] += p
        self._d += p

    # getters
    def get_bp_type(self):
        return self._bp_type

    def get_copy(self):
        return Basepair.copy(self)

    def get_d(self):
        return self._d

    def get_name(self):
        return self._name

    def get_partner(self, r_uuid):
        """
        get the other basepairing partner of a residue will throw an error
        if the supplied residue is not contained within this basepair

        :param res: the residue that you want to get the partner of
        :type res: secondary_structure.Residue object
        """

        if   r_uuid == self._res1_uuid:
            return self._res2_uuid
        elif r_uuid == self._res2_uuid:
            return self._res1_uuid
        else:
            raise exceptions.BasepairException(
                "call partner with a residue not in basepair")

    def get_r(self):
        return self._r

    def get_res1_sugar(self):
        return self._sugars[0]

    def get_res2_sugar(self):
        return self._sugars[1]

    def get_res1_uuid(self):
        return self._res1_uuid

    def get_res2_uuid(self):
        return self._res2_uuid

    def get_res_uuids(self):
        return [self._res1_uuid, self._res2_uuid]

    def get_state(self):
        sugars = [ np.copy(self._sugars[0]), np.copy(self._sugars[1])]
        return motif_state.Basepair(self._res1_uuid, self._res2_uuid, np.copy(self._r),
                                    np.copy(self._d), sugars, self._name, self._bp_type,
                                    self._uuid)

    def get_str(self):
        s  = basic_io.point_to_str(self._d) + ";"
        s += basic_io.matrix_to_str(self._r) + ";"
        s += basic_io.points_to_str(self._sugars) + ";"
        s += self._name + ";" + str(self._x3dna_bp_type) + ";" + str(self._bp_type) + ";"
        return s

    def get_sugars(self):
        return self._sugars

    def get_transforming_r_and_t(self, r, t, sugars):
        """
        get a rotation matrix and translation that describes the tranformation
        between the rotation, translation to THIS BasepairState.

        :param r: Another orientation matrix from another basepair
        :param t: The origin of another basepair
        :param sugars: the c1' coords of another basepair

        :type r: np.array
        :type t: np.array
        :type sugars: list of two np.arrays

        :return: rotation and translation that defines transformation betwen
            both states
        """

        r1 = self.r
        r2 = r
        r_trans = util.unitarize(r1.T.dot(r2))
        t_trans = -t

        new_sugars_2 = np.dot(sugars, r_trans.T) + t_trans + self.d

        if sugars is not None:
            diff = (((self.sugars[0] - new_sugars_2[0]) +
                     (self.sugars[1] - new_sugars_2[1]))/2)
        else:
            diff = 0
        return r_trans, t_trans+diff

    def get_transforming_r_and_t_w_state(self, state):
        """
        wrapper for get_transforming_r_and_t using another basepair state
        instead of specifying each component explicitly.

        :param state: the basepair state you would like to get a transformation
            to align to the basepair state calling this function
        :type state: BasepairState

        """

        return self.get_transforming_r_and_t(state._r, state._d, state._sugars)

    def get_transformed_state(self, r, t):
        """
        get new orientation, origin and sugar coordinates after transforming
        with suplied rotation and translation.

        :param r: supplied rotation matrix
        :param t: supplied translation

        :type r: np.array
        :type t: np.array

        :return: new orientation, origin and sugar coorindates of this basepair
            state after transformation

        """

        r_T = r.T

        new_r = util.unitarize(np.dot(self._r, r_T))
        new_sugars = np.dot(self._sugars, r_T) + t
        new_origin = np.dot(self._d, r_T) + t

        return new_r, new_origin, new_sugars

    def get_x3dna_bp_type(self):
        return self._x3dna_bp_type

    # private
    def _rot_diff(self, bp):
        r_diff = util.matrix_distance(self.r, bp.r)
        bp.flip()
        r_diff_2 = util.matrix_distance(self.r, bp.r)
        bp.flip()
        if r_diff > r_diff_2:
            r_diff = r_diff_2
        return r_diff


class RNAStructure(primitives.RNAStructure):
    """
    Complete container for representing a RNA. Contains both the 3D structure
    information but also includes basepair objects to represent the pairs
    between residues. RNAStructure is rarely called directly but serves as an
    abstract class for both Motif and Pose so for example use please see those
    classes.

    :param struct: The 3D coordinate information
    :param basepairs: basepairing information for each residue pair
    :param ends: the basepairs at the end of chains. These define connection
        points to other RNAStructures
    :param name: name of RNAStructure
    :param path: location of where RNAStructure where 3D information was loaded
        from either the pdb or directory.
    :param mtype: The type of a motif, only needed for Motif.motif objects
    :param score: The score generated by motif_scorer.MotifScorer

    :type struct: structure.Structure
    :type basepairs: list of basepair.Basepairs
    :type ends: list of basepair.Basepairs
    :type name: str
    :type path: str
    :type mtype: motif_type
    :type score: float

    :attributes:

    `structure`: structure.Structure
        structure containing residue and chain information for
        The 3D coordinate information
    `basepairs`: list of basepair.Basepairs
        Basepairs between residues
    `ends`: list of baspir.Basepairs
        Basepair ends where RNA structures can be connected
    `mtype`: motif_type
        The type of a motif, only needed for Motif.motif objects
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
        in the perspective of a given basepair end
    `beads` : list of residue.Bead objects
        keeps the beads required for steric clash calulations

    :examples:

    ..  code-block:: python

        # load test structure
        >>> from rnamake.unittests import instances
        >>> r_struct = instances.rna_structure()

    """

    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_score",
        "_protein_beads",
        "_end_ids",
        "_block_end_add",
        "_dot_bracket"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, name, dot_bracket,
                 block_end_add=-1, protein_beads=None):
        self._structure       = structure
        self._basepairs       = basepairs
        self._ends            = ends
        self._end_ids         = end_ids
        self._name            = name
        self._block_end_add   = block_end_add
        self._dot_bracket     = dot_bracket
        self._protein_beads   = protein_beads

        if self._protein_beads is None:
            self._protein_beads = []

        if len(self._end_ids) != len(self._ends):
            raise exceptions.RNAStructureException(
                "RNAStructures must have the same number of ends as end_ids "
                "has %d ends and % end_ids" % (len(self._ends), len(self._end_ids)))

        db_length = 0
        for e in self._dot_bracket:
            if e != "&":
                db_length += 1

        if db_length != self.get_num_res():
            raise exceptions.RNAStructureException(
                "RNAStructures must have dot_bracket symbol for each residue "
                "has %d residues and %d db symbols" % (self.get_num_res(), db_length))


        for c in self._structure.get_chains():
            for r in c:
                found = 0
                for i, end in enumerate(ends):
                    if (end.get_res1_uuid() == r.get_uuid() or
                        end.get_res2_uuid() == r.get_uuid()) and \
                        i != self._block_end_add:
                        found = 1
                        break
                if not found:
                    r.build_beads()

    def __iter__(self):
        return self._structure.__iter__()

    @classmethod
    def from_str(cls, s, rts):
        spl = s.split("&")
        name = spl[0]
        block_end_add = int(spl[1])
        struc = Structure.from_str(spl[2], rts)
        bp_strs = spl[3].split("@")
        bps = []
        for bp_str in bp_strs[:-1]:
            bps.append(bp_from_str(struc, bp_str))
        end_strs = spl[4].split("@")
        ends = []
        for end_str in end_strs[:-1]:
            ends.append(bp_from_str(struc, end_str))
        end_ids = spl[5].split()
        bead_strs = spl[6].split(";")
        protein_beads = []
        for bead_str in bead_strs[:-1]:
            protein_beads.append(bead.Bead.from_str(bead_str))

        dot_bracket = ""
        for i in range(7, len(spl) - 1):
            dot_bracket += spl[i]
            if i != len(spl) - 2:
                dot_bracket += "&"

        return cls(struc, bps, ends, end_ids, name, dot_bracket, block_end_add,
                   protein_beads)

    @classmethod
    def copy(cls, rs, new_uuid=0):
        s = Structure.copy(rs._structure, new_uuid)
        basepairs = []
        ends = []
        protein_beads = [ b.get_copy() for b in rs._protein_beads ]

        for bp in rs._basepairs:
            if new_uuid:
                bp_res = rs.get_bp_res(bp)
                res1 = s.get_residue(num=bp_res[0].get_num(),
                                     chain_id=bp_res[0].get_chain_id(),
                                     i_code=bp_res[0].get_i_code())
                res2 = s.get_residue(num=bp_res[1].get_num(),
                                     chain_id=bp_res[1].get_chain_id(),
                                     i_code=bp_res[1].get_i_code())
                bp_new = Basepair.copy_with_new_uuids(bp, res1.get_uuid(),
                                                      res2.get_uuid())
                basepairs.append(bp_new)
            else:
                basepairs.append(Basepair.copy(bp))

        for end in rs._ends:
            if new_uuid:
                bp_res = rs.get_bp_res(end)
                res1 = s.get_residue(num=bp_res[0].get_num(),
                                     chain_id=bp_res[0].get_chain_id(),
                                     i_code=bp_res[0].get_i_code())
                res2 = s.get_residue(num=bp_res[1].get_num(),
                                     chain_id=bp_res[1].get_chain_id(),
                                     i_code=bp_res[1].get_i_code())
                bp = Basepair.copy_with_new_uuids(end, res1.get_uuid(),
                                                  res2.get_uuid())
                ends.append(bp)
            else:
                ends.append(Basepair.copy(end))

        return cls(s, basepairs, ends, rs._end_ids, rs._name, rs._dot_bracket,
                   rs._block_end_add, protein_beads)

    def move(self, p):
        self._structure.move(p)
        for bp in self._basepairs:
            bp.move(p)
        for end in self._ends:
            end.move(p)

    def transform(self, t):
        self._structure.transform(t)
        for bp in self._basepairs:
            bp.transform(t)
        for end in self._ends:
            end.transform(t)

    def steric_clash(self, rs, clash_radius=settings.CLASH_RADIUS):
        for r1 in self:
            for b1 in r1.iter_beads():
                for r2 in rs:
                    for b2 in r2.iter_beads():
                        if b1.btype == bead.BeadType.PHOS or \
                           b2.btype == bead.BeadType.PHOS:
                            continue
                        dist = util.distance(b1.center, b2.center)
                        if dist < clash_radius:
                            return 1
        return 0

    # getters
    def get_block_end_add(self):
        return self._block_end_add

    def get_copy(self, new_uuid=0):
        return RNAStructure.copy(self, new_uuid)

    def get_dot_bracket(self):
        """
        secondary structure for rna_structure
         """

        return self._dot_bracket

    def get_pdb_str(self, renumber=-1, close_chain=0):
        """
        creates a PDB string formatted verision of this Structure object.

        :param renumber: what should the first residue be numbered. -1 is
            to NOT renumber, Default=-1.
        :param close_chain: fixes the phosphate backbone, takes a while, so
            default is 0 not to run

        :type renumber: int
        :type close_chain: int

        :return: str
        """
        if close_chain:
            for c in self.iter_chains():
                chain_closure.close_chain(c)

        return self._structure.get_pdb_str(renumber)

    def get_pdb(self, fname="motif.pdb", renumber=-1, close_chain=0):
        """
        write structure to pdb file

        :param fname: name of the file of the pdb file you want to write to
        :param renumber: what should the first residue be numbered. -1 is
            to NOT renumber, Default=-1.
        :param close_chain: fixes the phosphate backbone, takes a while, so
            default is 0 not to run

        :type fname: str
        :type renumber: int
        :type close_chain: int

        :return: None

        """
        #if close_chain:
        #    for c in self._chains():
        #        chain_closure.close_chain(c)

        return self._structure.get_pdb(fname, renumber)

    def get_sequence(self):
        """
        wrapper for :func:`rnamake.secondary_structure.Structure.sequence`
        """

        seq = ""
        for i, c in enumerate(self._structure.get_chains()):
            if i != 0:
                seq += "&"
            for r in c:
                seq += r.get_name()
        return seq

    def get_secondary_structure(self):
        db_chains = self._dot_bracket.split("&")
        res = []
        for i, c in enumerate(self._structure.get_chains()):
            for j, r in enumerate(c):
                s_r = secondary_structure.Residue(
                    r.get_name(), db_chains[i][j], r.get_num(),
                    r.get_chain_id(), r.get_i_code(), r.get_uuid())
                res.append(s_r)

        s = secondary_structure.Structure(res, self._structure._chain_cuts)
        bps = []
        for bp in self._basepairs:
            s_bp = secondary_structure.Basepair(
                        bp.get_res1_uuid(), bp.get_res2_uuid(),
                        bp.get_name(), bp.get_uuid())
            bps.append(s_bp)
        ends = []
        for end in self._ends:
            s_end = secondary_structure.Basepair(
                        end.get_res1_uuid(), end.get_res2_uuid(),
                        end.get_name(), end.get_uuid())
            ends.append(s_end)

        return secondary_structure.RNAStructure(s, bps, ends, self._end_ids[::],
                                                self._name)

    def get_str(self):
        """
        stringifies rna structure object
        """
        s  = self._name + "&" + str(self._block_end_add) + "&"
        s += self._structure.get_str() + "&"
        for bp in self._basepairs:
            res1, res2 = self.get_bp_res(bp)
            s += bp.get_str() + ";"
            s += str(res1.get_num()) + "|" + res1.get_chain_id() + "|" + \
                 res1.get_i_code() + ";"
            s += str(res2.get_num()) + "|" + res2.get_chain_id() + "|" + \
                 res2.get_i_code() + "@"
        s += "&"
        for end in self._ends:
            res1, res2 = self.get_bp_res(end)
            s += end.get_str() + ";"
            s += str(res1.get_num()) + "|" + res1.get_chain_id() + "|" + \
                 res1.get_i_code() + ";"
            s += str(res2.get_num()) + "|" + res2.get_chain_id() + "|" + \
                 res2.get_i_code() + "@"
        s += "&"
        for end_id in self._end_ids:
            s += end_id + " "
        s += "&"
        s += basic_io.beads_to_str(self._protein_beads)
        s += "&"
        s += self._dot_bracket
        s += "&"
        return s



class Motif(object):
    """
    The basic unit of this project stores the 3D coordinates of a RNA Motif
    as well as the 3DNA parameters such as reference frame and origin for
    each basepair

    :param mdir: the path to a motif directory that contains required files
    :type mdir: str

    :param pdb: the path to a pdb file to create this motif object, will also
        creates ref_frames.dat file and dssr out file in current directory
    :type pdb: str

    :param mtype: the enum motif type that this motif is, default UNKNOWN
    :type mtype: motif_type enum

    .. code-block:: python
        #creation from motif dir (recommended)

        #creation from a pdb, generates x3dna files at runtime
        >>> Motif(pdb="test.pdb")
        <Motif(name='test', ends='0')>

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
        that are at the end of Motif this is critical to assemble motifs
        together its not necessarily the first and last Basepairs
    `structure` : Structure object
        holds 3D coordinate data
    `name` : str
        the name of the directory without entire path
    """

    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_score",
        "_protein_beads",
        "_end_ids",
        "_block_end_add",
        "_dot_bracket",
        "_uuid",
        "_mtype",
        "_score"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, name, mtype, score,
                 dot_bracket, block_end_add=0, protein_beads=None, m_uuid=None):

        super(self.__class__, self).__init__(structure, basepairs, ends,
                                             end_ids, name, dot_bracket,
                                             block_end_add, protein_beads)

        self._mtype = mtype
        self._score = score
        self._uuid = m_uuid

        if self._uuid is None:
            self._uuid = uuid.uuid1()

    @classmethod
    def from_str(cls, s, rts):
        spl = s.split("&")
        name = spl[0]
        score = float(spl[1])
        block_end_add = int(spl[2])
        mtype = int(spl[3])
        struc = structure.Structure.from_str(spl[4], rts)
        bp_strs = spl[5].split("@")
        bps = []
        for bp_str in bp_strs[:-1]:
            bps.append(rna_structure.bp_from_str(struc, bp_str))
        ends = []
        end_strs = spl[6].split("@")
        for end_str in end_strs[:-1]:
            ends.append(rna_structure.bp_from_str(struc, end_str))
        end_ids = spl[7].split()
        bead_strs = spl[8].split(";")
        protein_beads = []
        for bead_str in bead_strs[:-1]:
            protein_beads.append(bead.Bead.from_str(bead_str))

        dot_bracket = ""
        for i in range(9, len(spl)-1):
            dot_bracket += spl[i]
            if i != len(spl)-2:
                dot_bracket += "&"

        return cls(struc, bps, ends, end_ids, name, mtype, score,
                   dot_bracket, block_end_add, protein_beads)

    @classmethod
    def copy(cls, m, new_uuid=0):
        s = structure.Structure.copy(m._structure, new_uuid)
        basepairs = []
        ends = []
        protein_beads = [ residue.Bead.copy(b) for b in m._protein_beads ]
        m_uuid = m._uuid

        for bp in m._basepairs:
            if new_uuid:
                bp_res = m.get_bp_res(bp)
                r_pos_1 = m.get_res_index(bp_res[0])
                r_pos_2 = m.get_res_index(bp_res[1])
                res1 = s.get_residue(index=r_pos_1)
                res2 = s.get_residue(index=r_pos_2)
                new_bp = basepair.Basepair.copy_with_new_uuids(bp, res1.uuid, res2.uuid)
                basepairs.append(new_bp)
            else:
                basepairs.append(basepair.Basepair.copy(bp))

        for end in m._ends:
            if new_uuid:
                bp_res = m.get_bp_res(end)
                r_pos_1 = m.get_res_index(bp_res[0])
                r_pos_2 = m.get_res_index(bp_res[1])
                res1 = s.get_residue(index=r_pos_1)
                res2 = s.get_residue(index=r_pos_2)
                new_end = basepair.Basepair.copy_with_new_uuids(end, res1.uuid, res2.uuid)
                ends.append(new_end)
            else:
                ends.append(basepair.Basepair.copy(end))

        if new_uuid:
            m_uuid = uuid.uuid1()

        return cls(s, basepairs, ends, m._end_ids, m._name, m._mtype, m._score,
                   m._dot_bracket, m._block_end_add, protein_beads, m_uuid)

    @classmethod
    def altered_copy(cls, m, name=None, mtype=None):
        if name is None:
            name = m._name
        if mtype is None:
            mtype = m._mtype

        s = structure.Structure.copy(m._structure, new_uuid)
        basepairs = []
        ends = []
        protein_beads = [ residue.Bead.copy(b) for b in m._protein_beads ]

        for bp in m._basepairs:
            if new_uuid:
                bp_res = m.get_bp_res(bp)
                res1 = s.get_residue(num=bp_res[0].num, chain_id=bp_res[0].chain_id,
                                     i_code=bp_res[0].i_code)
                res2 = s.get_residue(num=bp_res[1].num, chain_id=bp_res[1].chain_id,
                                     i_code=bp_res[1].i_code)
                bp = basepair.Basepair.copy_with_new_uuids(bp, res1.uuid, res2.uuid)
                basepairs.append(bp)
            else:
                basepairs.append(basepair.Basepair.copy(bp))

        for end in m._ends:
            i = m._basepairs.index(end)
            ends.append(basepairs[i])

        return cls(s, basepairs, ends, m._end_ids, name, mtype, m._score,
                   m._block_end_add, m._dot_bracket, protein_beads, m._uuid)

    def __repr__(self):
        """
        is called when motif is printed
        """
        return "<Motif(\n\tstructure='%s', \n\tends='%s')>" % (
        self._structure,len(self._ends))

    def to_str(self):
        """
        stringifies motif object
        """
        s  = self._name + "&" + str(self._score) + "&"
        s += str(self._block_end_add) + "&"
        s += str(self._mtype) + "&" + self._structure.to_str() + "&"
        for bp in self._basepairs:
            res1, res2 = self.get_bp_res(bp)
            s += bp.to_str() + ";"
            s += str(res1.num) + "|" + res1.chain_id + "|" + res1.i_code + ";"
            s += str(res2.num) + "|" + res2.chain_id + "|" + res2.i_code + "@"
        s += "&"
        for end in self._ends:
            res1, res2 = self.get_bp_res(end)
            s += end.to_str() + ";"
            s += str(res1.num) + "|" + res1.chain_id + "|" + res1.i_code + ";"
            s += str(res2.num) + "|" + res2.chain_id + "|" + res2.i_code + "@"
        s += "&"
        for end_id in self._end_ids:
            s += end_id + " "
        s += "&"
        s += basic_io.beads_to_str(self._protein_beads)
        s += "&"
        s += self._dot_bracket
        s += "&"
        return s

    def get_secondary_structure(self):
        db_chains = self._dot_bracket.split("&")
        res = []
        for i, c in enumerate(self._structure.get_chains()):
            for j, r in enumerate(c):
                s_r = secondary_structure.Residue(r.name, db_chains[i][j], r.num,
                                                  r.chain_id, r.i_code, r.uuid)
                res.append(s_r)
        s = secondary_structure.Structure(res, self._structure._chain_cuts)
        bps = []
        for bp in self._basepairs:
            s_bp = secondary_structure.Basepair(bp.res1_uuid, bp.res2_uuid,
                                                bp.name, bp.uuid)
            bps.append(s_bp)
        ends = []
        for end in self._ends:
            s_end = secondary_structure.Basepair(end.res1_uuid, end.res2_uuid,
                                                 end.name, end.uuid)
            ends.append(s_end)

        return secondary_structure.Motif(s, bps, ends, self._end_ids[::],
                                         self._mtype, self._name, self._uuid)

    def get_state(self):
        basepairs = []
        for bp in self._basepairs:
            sugars = [ np.copy(bp.sugars[0]), np.copy(bp.sugars[1])]
            bp_state = motif_state.Basepair(bp.res1_uuid, bp.res2_uuid, np.copy(bp.r),
                                            np.copy(bp.d), sugars, bp.name, bp.x3dna_bp_type,
                                            bp.bp_type, bp.uuid)
            basepairs.append(bp_state)

        ends = []
        for end in self._ends:
            sugars = [np.copy(end.sugars[0]), np.copy(end.sugars[1])]
            bp_state = motif_state.Basepair(end.res1_uuid, end.res2_uuid, np.copy(end.r),
                                            np.copy(end.d), sugars, end.name, end.x3dna_bp_type,
                                            end.bp_type, end.uuid)
            ends.append(bp_state)

        # keep only watson and crick bps
        #bps = []
        #for bp in basepairs:
        #    if bp.bp_type == "cW-W":
        #        bps.append(bp)
        ms = motif_state.Motif(self._structure.get_state(), basepairs, ends, self._end_ids,
                               self._name, self._mtype, self._score, self._dot_bracket,
                               self._block_end_add, self._uuid)
        return ms

    @property
    def score(self):
        return self._score

    @property
    def mtype(self):
        return self._mtype

    @property
    def uuid(self):
        return self._uuid


class MotifAligner(primitives.Aligner):
    __slots__ = [
        '_r',
        '_d',
        '_t',
        '_m_end',
        '_bp_pos_diff',
        '_dist_1',
        '_dist_2',
        '_sugar_diff_1',
        '_sugar_diff_2']

    def __init__(self):
        super(self.__class__, self).__init__()
        self._bp_pos_diff = 0
        self._m_end = None
        self._dist_1 = 0
        self._dist_2 = 0
        self._sugar_diff_1 = 0
        self._sugar_diff_2 = 0

    def get_aligned(self, ref_bp, m):
        m_copy = Motif.copy(m)
        self.align(ref_bp, m_copy)
        return m_copy

    def align(self, ref_bp, m):
        self._m_end = m.get_end(0)

        # calculate rotation before ref_bp (where we are going) and current
        # base pair (where we are)
        self._r = util.unitarize(ref_bp.r.T.dot(self._m_end.r))
        self._d = -self._m_end.d
        self._t = transform.Transform(self._r, self._d)
        m.transform(self._t)
        self._bp_pos_diff = ref_bp.d - self._m_end.d
        m.move(self._bp_pos_diff)

        # alignment is by center of basepair, it can be slightly improved by
        # aligning the c1' sugars
        self._dist_1 = util.distance(self._m_end.res1_sugar, ref_bp.res1_sugar)
        self._dist_2 = util.distance(self._m_end.res2_sugar, ref_bp.res1_sugar)

        if self._dist_1 < self._dist_2:
            self._sugar_diff_1 = ref_bp.res1_sugar - self._m_end.res1_sugar
            self._sugar_diff_2 = ref_bp.res2_sugar - self._m_end.res2_sugar
        else:
            self._sugar_diff_1 = ref_bp.res1_sugar - self._m_end.res2_sugar
            self._sugar_diff_2 = ref_bp.res2_sugar - self._m_end.res1_sugar

        if self._dist_1 < 5 or self._dist_2 < 5:
            m.move((self._sugar_diff_1 + self._sugar_diff_2) / 2)


# DEPRECATED
def align_motif(ref_bp_state, motif_end, motif):
    """
    This is the workhorse of the entire suite. Aligns one end of a motif to
    the reference frame and origin of a Basepair object.

    :param ref_bp: the base pair that the motif end is going to align too
    :param motif_end: the motif end basepair to overly with the ref_bp
    :param motif: the motif object that you want to align

    :type ref_bp: Basepair object
    :type motif_end: Basepair object
    :type motif: Motif object
    """

    r1 , r2 = ref_bp_state.r , motif_end.r
    r = util.unitarize(r1.T.dot(r2))
    trans = -motif_end.d
    t = transform.Transform(r, trans)
    motif.transform(t)
    bp_pos_diff = ref_bp_state.d - motif_end.d
    motif.move(bp_pos_diff)

    #alignment is by center of basepair, it can be slightly improved by
    #aligning the c1' sugars
    res1_coord, res2_coord = motif_end.sugars
    ref_res1_coord, ref_res2_coord = ref_bp_state.sugars

    dist1 = util.distance(res1_coord, ref_res1_coord)
    dist2 = util.distance(res2_coord, ref_res1_coord)

    if dist1 < dist2:
        sugar_diff_1 = ref_res1_coord - res1_coord
        sugar_diff_2 = ref_res2_coord - res2_coord
    else:
        sugar_diff_1 = ref_res1_coord - res2_coord
        sugar_diff_2 = ref_res2_coord - res1_coord

    if dist1 < 5 or dist2 < 5:
        motif.move( (sugar_diff_1 + sugar_diff_2) / 2 )

# DEPRECATED
def get_aligned_motif(ref_bp, motif_end, m, sterics=1):

    motif_end_index = m.get_end_index(motif_end.name)
    m_copy = Motif.copy(m)
    motif_end = m_copy.get_end(motif_end_index)
    align_motif(ref_bp, motif_end, m_copy)

    return m_copy


def are_residues_connected(r1, r2, cutoff=3.0):
        """
        Determine if another residue is connected to this residue, returns 0
        if res is not connected to self, returns 1 if connection is going
        from 5' to 3' and returns -1 if connection is going from 3' to 5'

        :param res: another residue
        :param cutoff: distance to be considered connected, default: 3 Angstroms

        :type res: Residue
        :type cutoff: int

        :rtype: int

        """
        # 5' to 3'
        if r1.has_atom("O3'") and r2.has_atom("P"):
            o3_atom = r1.get_atom("O3'")
            p_atom = r2.get_atom("P")
            if util.distance(o3_atom.get_coords(), p_atom.get_coords()) < cutoff:
                return 1

        # 3' to 5'
        if r1.has_atom("P") and r2.has_atom("O3'"):
            p_atom = r1.get_atom("P")
            o3_atom = r2.get_atom("O3'")
            if util.distance(o3_atom.get_coords(), p_atom.get_coords()) < cutoff:
                return -1

        return 0


def connect_residues_into_chains(residues):
        """
        takes all residues and puts into the correct order in chains checking
        for physical connection between O5' and P atoms between residues

        :param residues: residue objects that belong in this structure
        :type residues: List of Residue objects

        :return: List of Chain objects
        """

        chains = []
        # sort residues so check residues for connection quicker as the next on
        # in the array will be closest to it by number
        residues.sort(key=lambda x: x.get_num())

        while True:
            current = None
            # find next 5' end, all chains go from 5' to 3'
            for i, r in enumerate(residues):
                five_prime_end = 1
                for j, r2 in enumerate(residues):
                    if are_residues_connected(r, r2) == -1:
                        five_prime_end = 0
                        break
                if five_prime_end:
                    current = r
                    break
            if not current:
                break
            residues.remove(current)
            current_chain_res = []
            # extend chain until 3' end
            while current is not None:
                current_chain_res.append(current)
                found = 0
                for r in residues:
                    if are_residues_connected(current, r) == 1:
                        current = r
                        found = 1
                        break
                if found:
                    residues.remove(current)
                else:
                    # no more residues to add, make chain object
                    chains.append(Chain(current_chain_res))
                    current = None

        return chains


def _calc_center(res):
    center = np.array([0.0,0.0,0.0])
    count = 0
    for r in res:
        for a in r:
            if a is None:
                continue
            center += a.get_coords()
            count += 1
    center /= count
    return center


def bp_from_str(struc, s):
    bp_spl = s.split(";")
    r2_info = bp_spl.pop().split("|")
    r1_info = bp_spl.pop().split("|")
    bp_str = ";".join(bp_spl)
    res1 = struc.get_residue(int(r1_info[0]), r1_info[1], r1_info[2])
    res2 = struc.get_residue(int(r2_info[0]), r2_info[1], r2_info[2])
    return Basepair.from_str(bp_str, res1.get_uuid(), res2.get_uuid())


def basepairs_from_x3dna(path, s):
    """
    gets x3dna data on basepairing information and then interwines it
    with the structural information stored in structure for simpler
    retreival of data

    :param path: path to the pdb file
    :param structure: the structure with the same residues that will appear
        in the x3dna output

    :type path: str
    :type structure: structure.Structure

    :return: gets all the basepairs that x3dna finds and returns them as
        basepair.Basepair
    :rtype: list of basepair.Basepair objects

    """
    x3dna_parser = x3dna.X3dna()
    x_basepairs = x3dna_parser.get_basepairs(path)
    basepairs = []
    for xbp in x_basepairs:
        res1 = s.get_residue(num=xbp.res1.num,
                             chain_id=xbp.res1.chain_id,
                             i_code=xbp.res1.i_code)

        res2 = s.get_residue(num=xbp.res2.num,
                             chain_id=xbp.res2.chain_id,
                             i_code=xbp.res2.i_code)

        if res1 is None or res2 is None:
            not_found = 0
            if res1 is None:
                not_found = xbp.res1
            else:
                not_found = xbp.res2
            user_warnings.RNAStructureWarning(
                "cannot find residues in basepair: " + str(not_found) + " "
                "this residue should NOT be a normal nucleotide etc: A,G,C,U\n")
            continue

        try:
            res1.get_atom("C1'")
        except exceptions.ResidueException:
            user_warnings.RNAStructureWarning(
                str(res1) + " has no C1' residue cannot have it in a basepair "
                "is required for alignment\n")
            continue

        try:
            res2.get_atom("C1'")
        except exceptions.ResidueException:
            user_warnings.RNAStructureWarning(
                str(res2) + " has no C1' residue cannot have it in a basepair "
                "is required for alignment\n")
            continue


        center = _calc_center([res1, res2])
        name = primitives.calc_bp_name([res1, res2])
        bp_type = primitives.BasepairType.NC

        bp_str = res1.get_name()+res2.get_name()
        wc = "GC,CG,AU,UA".split(",")

        if bp_str in wc and xbp.bp_type == x3dna.X3dnaBPType.cWUW:
            bp_type = primitives.BasepairType.WC
        elif bp_str == "GU" or bp_str == "UG" and xbp.bp_type == x3dna.X3dnaBPType.cWUW:
            bp_type = primitives.BasepairType.GU

        bp = Basepair(res1.get_uuid(), res2.get_uuid(), xbp.r, center,
                      [res1.get_coords("C1'"), res2.get_coords("C1'")],
                      name, bp_type=bp_type, x3dna_bp_type=xbp.bp_type)
        basepairs.append(bp)

    """if os.path.isfile("ref_frames.dat"):
        os.remove("ref_frames.dat")

    if os.path.isfile(name + "_dssr.out"):
        os.remove(name + "_dssr.out")"""

    return basepairs


def get_chain_end_map(chains, end):
    """
    builds :class:`ChainEndPairMap` instances from a target end and a pool
    of chains from a structure.

    :param chains: all chains that could contain residues in the target end
    :param end:

    :type chains: list of chain.Chain objects
    :type end: basepair.Basepair

    :return: a ChainEndPairMap defining the 5' and 3' chains of a given end
    :rtype: ChainEndPairMap
    """

    chain_map = ChainEndPairMap()

    for c in chains:
        for r in end.residues():
            if c.first().uuid == r.uuid and chain_map.p5_chain is None:
                chain_map.p5_chain = c
            elif c.first().uuid == r.uuid:
                raise exceptions.RNAStructureException(
                    "cannot build chain map two residues are assigned to 5' "
                    "chain")

            if c.last().uuid == r.uuid and chain_map.p3_chain is None:
                chain_map.p3_chain = c
            elif c.last().uuid == r.uuid:
                raise exceptions.RNAStructureException(
                    "cannot build chain map two residues are assigned to 3' "
                    "chain")

    if chain_map.p5_chain is None or chain_map.p3_chain is None:
        raise exceptions.RNAStructureException(
            "did not build map properly, both chains are not found")

    return chain_map


def calc_center(res):
    center = np.array([0.0,0.0,0.0])
    count = 0
    for r in res:
        for a in r:
            if a is None:
                continue
            center += a.get_coords()
            count += 1
    center /= count
    return center


# TODO should probably move this somewhere else?
class ChainEndPairMap(object):
    """
    A simple class container for expressing the relationship between chain
    objects and end objects. An end will always have two residues one at the
    end of a 5' end of a chain and the other at the 3' end of a chain. This
    class organizes the two chains which may be different chains to match the
    two residues in an end. This is primarily used when mergering chains
    together to get the final sequence of a construct.

    :param chain1: the chain that res1 is at the 5' end at in the basepair that this
        object is mapping for
    :param chain2: the chain that res1 is at the 3' end at in the basepair that this
        object is mapping for

    :type chain1: chain.Chain
    :type chain2: chain.Chain

    :attributes:

    `chain1`: chain.Chain
        the chain that res1 is at the 5' end at in the basepair that this
        object is mapping for
    `chain2`: chain.Chain
     the chain that res1 is at the 3' end at in the basepair that this
        object is mapping for

    """

    def __init__(self, chain1=None, chain2=None):
        self.p5_chain, self.p3_chain = chain1, chain2

    def is_hairpin(self):
        """
        A quick check to see if both chains are the same. This would mean
        the chain is a hairpin since both chain ends form a basepair

        :returns: whether both chains are the same
        :rtype: int

        """

        if self.p5_chain == self.p3_chain:
            return 1
        else:
            return 0

    def chains(self):
        """
        returns both chains into an list for easy iteration

        :returns: both chains in a list
        :rtype: list of chain.Chains
        """

        return [self.p5_chain, self.p3_chain]


