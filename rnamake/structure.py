import numpy as np

import pdb_parser
import chain
import util
import exceptions

class Structure(object):
    """Stores 3D structure information from a pdb file. Stores all chains,
    residues and atoms objects. Implementation is designed to be extremely
    lightweight and capable of performing fast transformations. to load a PDB
    formated file into a Structure object use
    :func:`structure_from_pdb`

    :param chains: chains that belong to this structure, optional
    :type chains: list of chain.Chain objects

    :attributes:
    `chains` : list of Chain objects that belong to the current structure

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

    def __init__(self, chains=None):
        self.chains = chains
        if self.chains is None:
            self.chains = []
        self.name = "N/A"

    def __repr__(self):
        return """<Structure(name: %s, #chains: %s, #residues: %s, #atoms: %s)>""" %\
               (self.name, len(self.chains), len(self.residues()), len(self.atoms()))

    def get_beads(self, excluded_res=None):
        """
        generates 3-bead model residue beads for all residues in current
        structure.

        :param excluded_res: List of residue objects whose beads are not to be
            included. This is generally end residues that would instantly clash
            with residues they are being overlayed onto when performing motif
            aligning
        :type excluded_res: List of Residue objects

        :return: List of Bead objects

        :examples:

        .. code-block:: python

            # load structure from pdb formatted file
            >>> import rnamake.unittests.files
            >>> s = structure_from_pdb(rnamake.unittests.files.P4P6_PDB_PATH)

            >>> len(s.get_beads())
            470

            >>> s.get_beads()[0]
            <Bead(btype='SUGAR', center='-24.027 -48.5001111111 86.368')>

            # exclude beads from first residue. This can be useful if you only
            # need sterics from a part of the structure
            >>> len(s.get_beads(excluded_res=[s.residues()[0]]))
            468
        """

        if excluded_res is None:
            excluded_res = []

        beads = []
        for r in self.residues():
            if r in excluded_res:
                continue
            beads.extend(r.get_beads())
        return beads

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

              # load structure from pdb formatted file
              >>> import rnamake.unittests.files
              >>> s = structure_from_pdb(rnamake.unittests.files.P4P6_PDB_PATH)

              # get specific residue
              >>> r = s.get_residue(num=106)
              >>> print r
              <Residue('U106 chain A')>

              # get residue using its unique indentifer
              >>>s.get_residue(uuid=r.uuid)
              <Residue('U106 chain A')>
        """

        #nothing specified
        if num is None and chain_id is None and i_code is None and uuid is None:
            raise exceptions.StructureException("called get_residue wiht no "
                                                "arguments")

        found = []
        for c in self.chains:
            for r in c.residues:
                if uuid is not None and uuid != r.uuid:
                    continue
                if num is not None and num != r.num:
                    continue
                if i_code is not None and i_code != r.i_code:
                    continue
                if chain_id is not None and chain_id != r.chain_id:
                    continue
                found.append(r)

        if len(found) > 1:
            raise exceptions.StructureException(
                "found multiple residues in get_residue(), narrow " +
                "your search")

        # TODO maybe add warning here?
        if len(found) == 0:
            return None

        return found[0]

    def residues(self):
        """
        Concats all residue objects from all Chain objects intos a unified
        list to be able to easily iterate through.

        :return: List of Residue objects
        """
        residues = []
        for c in self.chains:
            residues.extend(c.residues)
        return residues

    def atoms(self):
        """
        Concats all Atom objects from all Residue objects intos a unified
        list to be able to easily iterate through.

        :return: List of Residue objects
        """

        atoms = []
        for r in self.residues():
            for a in r.atoms:
                if a is None:
                    continue
                atoms.append(a)
        return atoms

    def to_str(self):
        """
        Stringifes Structure object

        :return: str
        """

        s = ""
        for c in self.chains:
            s += c.to_str() + ":"
        return s

    def to_pdb_str(self, renumber=-1):
        """
        creates a PDB string formatted verision of this Structure object.

        :param renumber: what should the first residue be numbered. -1 is
            to NOT renumber, Default=-1.
        :return: int

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

        for i, c in enumerate(self.chains):
            c_str, acount = c.to_pdb_str(acount, 1, rnum, chain_id)
            if renumber != -1:
            #    chain_id = c_names[i+1]
                rnum += len(c.residues)
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

    def copy(self):
        """
        creates a deep copy of this structure

        :returns: copy of Structure object
        :rtype: Structure
        """
        chains = []
        for c in self.chains:
            cc = c.copy()
            chains.append(cc)

        return Structure(chains)

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



        r_T = t.rotation().T
        for a in self.atoms():
            a.coords = np.dot(a.coords, r_T) + t.translation()

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

        for a in self.atoms():
            a.coords += p


def structure_from_pdb(pdb_path):
    """
    Processes a PDB formatted into Structure object. Uses pdb_parser module
    to accomplish this.

    :param pdb_path: path to PDB formatted file
    :type pdb_path: str

    :return: Structure object
    :rtype: Structure
    """

    residues = pdb_parser.parse(pdb_path)
    chains = chain.connect_residues_into_chains(residues)
    s = Structure(chains)
    s.name = util.filename(pdb_path[:-4])
    return s
