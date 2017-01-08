import uuid
import numpy as np

from atom import Atom
from bead import BeadType, Bead
import motif_state
import primitives.residue
import residue_type, util, basic_io, exceptions


class Residue(primitives.residue.Residue):
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

    `atoms` : atom.Atom
        holds all atoms that belong to this residue object
    `name` : str
        name of residue, ex. ADE, GUA etc
    `num` : int
        residue num
    `rtype` : residue_type.ResidueType
        Information about residue type each nucleic acid has its own type
    `chain_id` : str
        chain indentification string, ex. 'A' or 'B'
    `i_code`: str
        residue insertion code
    `uuid` : uuid.uuid1
        the unique id to indentify residue

    :examples:

    .. code-block:: python

        # generating a new residue
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
    def copy(cls, r, new_uuid=0):
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

        return cls(atoms, r._rtype, r._name, r._num, r._chain_id,
                   r._i_code, r_uuid)

    def __iter__(self):
        return self._atoms.__iter__()

    def __repr__(self):
        return "<Residue('%s%d%s chain %s')>" % (
            self._name, self._num, self._i_code, self._chain_id)

    def iter_beads(self):
        return self._beads.__iter__()

    def connected_to(self, res, cutoff=3.0):
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
        if self.has_atom("O3'") and res.has_atom("P"):
            o3_atom = self.get_atom("O3'")
            p_atom = res.get_atom("P")
            if util.distance(o3_atom.coords, p_atom.coords) < cutoff:
                return 1

        # 3' to 5'
        if self.has_atom("P") and res.has_atom("O3'"):
            p_atom = self.get_atom("P")
            o3_atom = res.get_atom("O3'")
            if util.distance(o3_atom.coords, p_atom.coords) < cutoff:
                return -1

        return 0

        return copied_r

    def to_str(self):
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
                s += a.to_str() + ","
        return s

    def to_pdb_str(self, acount=1, return_acount=0, rnum=-1, chain_id=""):
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
                 ('ATOM', acount, a.name, '', self.short_name(), cid,
                  num, '', a.coords[0], a.coords[1], a.coords[2], 1.00,
                  0.00, '', '')
            acount += 1

        if return_acount:
            return s, acount
        else:
            return s

    def to_pdb(self, fname="residue.pdb"):
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

    def has_atom(self, atom_name=None, index=None):
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

    def center(self):
        return util.center(self._atoms)

    def transform(self, t):
        for a in self._atoms:
            if a is not None:
                a.transform(t)
        for b in self._beads:
            b.transform(t)

    def move(self, p):
        for a in self._atoms:
            if a is not None:
                a.move(p)
        for b in self._beads:
            b.move(p)

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

    def get_coords(self, atom_name=None, index=None):
        return self.get_atom(atom_name, index).coords

    def short_name(self):
        """gets letter of residue, i.e. A or G etc

        :return: letter for residue
        :rtype: str
        """
        return self._rtype.short_name

    def num_beads(self):
        return len(self._beads)

    def get_state(self):
        beads = [Bead.copy(b) for b in self._beads]
        return motif_state.Residue(self._name, self._num, self._chain_id,
                                   self._i_code, beads, self._uuid)

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
                a = Atom(name_change[1], a.coords)
            if self._rtype.is_valid_atom(a.name):
                pos = self._rtype.atom_index(a.name)
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
        types = [BeadType.PHOS, BeadType.SUGAR, BeadType.BASE]
        for i, alist in enumerate([phos_atoms, sugar_atoms, base_atoms]):
            if len(alist) > 0:
                beads.append(Bead(util.center(alist), types[i]))

        return beads

    def build_beads(self):
        self._beads = self.__get_beads()
