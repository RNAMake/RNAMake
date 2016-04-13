import uuid
import numpy as np

import residue_type, util, basic_io

class BeadType(object):
    """
    BeadType is an ENUM type. This is to specify which center of atoms each bead
    represents.

    Phosphate (0):  P, OP1, OP2\n
    Sugar (1):  O5',C5',C4',O4',C3',O3',C1',C2',O2'\n
    Base  (2):  All remaining atoms
    """

    PHOS = 0   # P, OP1, OP2
    SUGAR = 1  # O5', C5', C4', O4', C3', O3', C1', C2', O2'
    BASE = 2   # All remaining


class Bead(object):
    """
    Bead class stores information related to keeping track of steric clashes
    between residues during building. They are never used outside the Residue class

    :param btype: type of the bead either Phos(Phosphate), Sugar or Base, of
        the atoms used generate the center
    :type btype: BeadType

    :param center: The geometric center of the group of atoms
    :type center: numpy array

    """

    __slots__ = ["center", "btype"]

    def __init__(self, center, btype):
        self.center, self.btype = center, btype

    def __repr__(self):
        center = basic_io.point_to_str(self.center)
        return "<Bead(btype='%s', center='%s')>" % (self.type_name(), center)

    def copy(self):
        """
        returns a deep copy of current bead object

        :returns: copy of bead object
        :rtype: Bead
        """
        return Bead(np.copy(self.center), self.btype)

    def type_name(self):
        """
        returns name of btype in string form

        :returns: name of btype type
        :rtype: str
        """
        if self.btype == 0:
            return "PHOSPHATE"
        if self.btype == 1:
            return "SUGAR"
        if self.btype == 2:
            return "BASE"


class Residue(object):
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
        >>>rts = residue_type.ResidueTypeSet()
        >>>rtype = rts.get_rtype_by_resname("ADE")
        >>>r = Residue(rtype, "ADE", 1, "A")
        >>>print r.name
        A

        #using test residue
        >>>import rnamake.unittests.instances
        >>>r = rnamake.unittests.instances.residue()
        >>>print r.name
        G

        >>>a = r.get_atom("C1'")
        >>>print a.coords
        [-23.806 -50.289  86.732]

        >>>r.get_beads()
        [<Bead(btype='SUGAR', center='-24.027 -48.5001111111 86.368')>, <Bead(btype='BASE', center='-21.2186363636 -52.048 85.1157272727')>]

        #a fast way of saving coordinate information to file
        >>>r.to_str()
        "GUA,G,103,A,,N,N,N,O5' -26.469 -47.756 84.669,C5' -25.05 -47.579 84.775,C4' -24.521 -48.156 86.068,O4' -24.861 -49.568 86.118,C3' -23.009 -48.119 86.281,O3' -22.548 -46.872 86.808,C1' -23.806 -50.289 86.732,C2' -22.812 -49.259 87.269,O2' -23.167 -48.903 88.592,N1 -19.538 -52.485 85.025,C2 -19.717 -51.643 86.097,N2 -18.624 -51.354 86.809,N3 -20.884 -51.124 86.445,C4 -21.881 -51.521 85.623,C5 -21.811 -52.356 84.527,C6 -20.546 -52.91 84.164,O6 -20.273 -53.677 83.228,N7 -23.063 -52.513 83.947,C8 -23.858 -51.786 84.686,N9 -23.21 -51.159 85.722,"

        #get PDB formmated coordinates back out
        >>>r.to_pdb_str()
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
        "rtype",
        "name",
        "num",
        "chain_id",
        "i_code",
        "uuid",
        "atoms"]

    def __init__(self, rtype, name, num, chain_id, i_code=""):
        self.rtype = rtype
        self.name = name
        self.num = num
        self.chain_id = chain_id
        self.i_code = i_code
        self.uuid = uuid.uuid1()
        # fix pickle bug
        if len(self.i_code) == 0:
            self.i_code = ""

        self.atoms = []

    def __repr__(self):
        return "<Residue('%s%d%s chain %s')>" % (
            self.name, self.num, self.i_code, self.chain_id)

    def short_name(self):
        return self.rtype.name[0]

    def setup_atoms(self, atoms):
        """
        put atoms in correct positon in internal atom list, also corrects some
        named atom names to their correct name

        :param atoms: list of atom objects that are to be part of this residue
        :type atoms: list of Atom objects
        """

        self.atoms = [None for x in self.rtype.atom_map.keys()]
        for a in atoms:
            name_change = self.rtype.get_correct_atom_name(a)
            if name_change is not None:
                a.name = name_change[1]
            if a.name in self.rtype.atom_map:
                pos = self.rtype.atom_map[a.name]
                self.atoms[pos] = a

    def get_atom(self, atom_name):
        """
        get atom object by its name

        :param atom_name: name
        :type atom_name: str

        Examples:

        .. code-block:: python

            >>>r = rnamake.unittests.instances.residue()
            >>>a = r.get_atom("C1'")
            >>>print a.coords
            [-23.806 -50.289  86.732]
        """

        try:
            index = self.rtype.atom_map[atom_name]
            return self.atoms[index]
        except KeyError:
            raise KeyError("cannot find atom")

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
        o3_atom = self.get_atom("O3'")
        p_atom = res.get_atom("P")
        if o3_atom and p_atom:
            if util.distance(o3_atom.coords, p_atom.coords) < cutoff:
                return 1

        # 3' to 5'
        p_atom = self.get_atom("P")
        o3_atom = res.get_atom("O3'")
        if o3_atom and p_atom:
            if util.distance(o3_atom.coords, p_atom.coords) < cutoff:
                return -1

        return 0

    def get_beads(self):
        """
        Generates steric beads required for checking for steric clashes between
        motifs. Each residues has three beads modeled after the typical three
        bead models used in coarse grain modeling. The three beads are:

        Phosphate:  P, OP1, OP2\n
        Sugar    :  O5',C5',C4',O4',C3',O3',C1',C2',O2'\n
        Base     :  All remaining atoms

        if there are for example no phosphate atoms only 2 beads will be returned.

        .. code-block:: python

            >>>import rnamake.unittests.instances
            >>>r = rnamake.unittests.instances.residue()
            >>>r.get_beads()
            [<Bead(btype='SUGAR', center='-24.027 -48.5001111111 86.368')>, <Bead(btype='BASE', center='-21.2186363636 -52.048 85.1157272727')>]

        """
        phos_atoms, sugar_atoms, base_atoms = [], [], []

        for i, a in enumerate(self.atoms):
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

    def copy(self):
        """
        performs a deep copy of Residue object

        :rtype: Residue

        :examples:

        .. code-block:: python

            >>>import rnamake.unittests.instances
            >>>r = rnamake.unittests.instances.residue()
            >>>r_copy = r.copy()
            >>>r_copy.name
            G

        """
        copied_r = Residue(self.rtype, self.name, self.num, self.chain_id, self.i_code)
        copied_r.atoms = [None for x in range(len(self.atoms))]
        for i, a in enumerate(self.atoms):
            if a is None:
                continue
            copied_r.atoms[i] = a.copy()

        copied_r.uuid = self.uuid
        return copied_r

    def new_uuid(self):
        """
        give residue a new unique indentifier code.
        There is probably no reason why you should call this unless writing a new,
        motif structure.
        """
        self.uuid = uuid.uuid1()

    def to_str(self):
        """
        stringifes residue object

        :returns: stringified residue object

        .. code-block:: python

            >>>import rnamake.unittests.instances
            >>>r = rnamake.unittests.instances.residue()
            >>>r.to_str()
            "GUA,G,103,A,,N,N,N,O5' -26.469 -47.756 84.669,C5' -25.05 -47.579 84.775,C4' -24.521 -48.156 86.068,O4' -24.861 -49.568 86.118,C3' -23.009 -48.119 86.281,O3' -22.548 -46.872 86.808,C1' -23.806 -50.289 86.732,C2' -22.812 -49.259 87.269,O2' -23.167 -48.903 88.592,N1 -19.538 -52.485 85.025,C2 -19.717 -51.643 86.097,N2 -18.624 -51.354 86.809,N3 -20.884 -51.124 86.445,C4 -21.881 -51.521 85.623,C5 -21.811 -52.356 84.527,C6 -20.546 -52.91 84.164,O6 -20.273 -53.677 83.228,N7 -23.063 -52.513 83.947,C8 -23.858 -51.786 84.686,N9 -23.21 -51.159 85.722,"
        """
        s = self.rtype.name + "," + self.name + "," + str(self.num) + "," + \
            self.chain_id + "," + self.i_code + ","
        for a in self.atoms:
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
            >>>r.to_pdb_str()
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

        num = self.num
        cid = self.chain_id
        if rnum != -1:
            num = rnum
        if chain_id != "":
            cid = chain_id

        s = ""
        for a in self.atoms:
            if a is None:
                continue
            s += basic_io.PDBLINE_GE100K % \
                 ('ATOM', acount, a.name, '', self.rtype.name[0], cid,
                  num, '', a.coords[0], a.coords[1], a.coords[2], 1.00,
                  0.00, '', '')
            acount += 1

        if return_acount:
            return s, acount
        else:
            return s

    def to_pdb(self, fname="residue.pdb"):
        f = open(fname, "w")
        s = self.to_pdb_str()
        f.write(s)
        f.close()



