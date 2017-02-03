import numpy as np

import pdb_parser
import util
import exceptions
import motif_state
import primitives.structure
from chain import Chain, connect_residues_into_chains
from residue import Residue

class Structure(primitives.structure.Structure):
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
               (len(self._chains), len(self._residues))

    def to_str(self):
        """
        Stringifes Structure object

        :return: str
        """

        s = ""
        for r in self._residues:
            s += r.to_str() + ";"
        s += " ".join([str(x) for x in self._chain_cuts]) + ";"
        return s

    def to_pdb_str(self, renumber=-1):
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
            c_str, acount = c.to_pdb_str(acount, 1, rnum, chain_id)
            if renumber != -1:
            #    chain_id = c_names[i+1]
                rnum += len(c)
            s += c_str
            s += "TER\n"
        return s

    def to_pdb(self, fname="structure.pdb", renumber=-1):
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

    def get_state(self):
        residues = [ r.get_state() for r in self._residues ]
        return motif_state.Structure(residues, self._chain_cuts)

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

def structure_from_pdb(pdb_path, rts):
    """
    Processes a PDB formatted into Structure object. Uses pdb_parser module
    to accomplish this.

    :param pdb_path: path to PDB formatted file
    :type pdb_path: str

    :return: Structure object
    :rtype: Structure
    """

    residues = pdb_parser.parse(pdb_path, rts=rts)

    chains = connect_residues_into_chains(residues)

    # TODO should probably clean this up at some point
    chain_cuts = []
    residues = []
    for i, c in enumerate(chains):
        for r in c:
            residues.append(r)
        if i < len(chains) - 1:
            chain_cuts.append(len(residues))
    chain_cuts.append(len(residues))

    s = Structure(residues, chain_cuts)
    return s
