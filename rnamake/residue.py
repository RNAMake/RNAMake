import uuid
import logging
from . import atom
from . import residue_type
from . import util

logging.basicConfig()
logger = logging.getLogger(__name__)

class BeadType(object):

    """
    BeadType is an ENUM type. This is to keep track what atoms are used to
    generate the bead center. There are three types Phos(Phosphate), Sugar and
    Base
    """
    Phos = 0
    Sugar = 1
    Base = 2


class Bead(object):

    """
    Bead class stores information related to keeping track of steric clashes
    between residues during building. They are never used outside the Residue class

    :param type: type of the bead either Phos(Phosphate), Sugar or Base, of
        the atoms used generate the center
    :btype btype: BeadType

    :param center: The geometric center of the group of atoms
    :type center: numpy array

    """

    __slots__ = ["center", "btype"]

    def __init__(self, center, btype):
        self.center, self.btype = center, btype


class Residue(object):

    """
    Store residue information from pdb file, stores all Atom objects that
    belong to residue. Implementation is designed to be extremely lightweight.

    :param rtype: residue type, stores information about residue
    :type rtype: ResidueType object

    :param name: residue name
    :type name: str

    :param num: residue num
    :type num: int

    :param chain_id: chain identification
    :type chain_id: str

    :param i_code: insertion code, optional argument if not supplied will be set to ""
    :type i_code: str
    """

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
        self.beads = []

    def __repr__(self):
        return "<Residue('%s%d%s chain %s')>" % (
            self.name, self.num, self.i_code, self.chain_id)

    def setup_atoms(self, atoms):
        """
        put atoms in correct positon in internal atom list, warns of extra atoms
        are included as well if atoms are missing

        :param atoms: list of atom objects that are to be part of this residue
        :type atoms: list of Atom objects
        """

        self.atoms = [None for x in self.rtype.atom_map.keys()]
        for a in atoms:
            #self.rtype.fix_common_alt_name(a)
            if a.name in self.rtype.atom_map:
                pos = self.rtype.atom_map[a.name]
                self.atoms[pos] = a
            else:
                logger.warning(a.name + " not included in " + repr(self))

        # Warning if an atom is missing
        for i,a in enumerate(self.atoms):
            if a is None:
                correct_name = None
                for name, pos in self.rtype.atom_map.iteritems():
                    if pos == i:
                        correct_name = name

                logger.warning(correct_name + " is undefined in " + repr(self))

    def get_atom(self, atom_name):
        """
        get atom object by its name

        :param atom_name: name
        :type atom_name: str
        .. code-block:: python
            >>>r.get_atom("P")
            <Atom(name='P', coords='0 1 2')>
        """
        try:
            index = self.rtype.atom_map[atom_name]
            return self.atoms[index]
        except KeyError:
            logger.critical("cannot find atom " + atom_name + "in " + repr(self))
            raise KeyError()
    def conntected_to(self,res,cutoff=3.0):
        """
		Determine if another residue is connected to this residue, returns 0
		if res is not connected to self, returns 1 if connection is going
		from 5' to 3' and returns -1 if connection is going from 3' to 5'

        :param res: another residue

        """
        pass

    def to_str(self):
        """
        stringifes residue object to string
        """
        s = self.rtype.name + "," + self.name + "," + str(self.num) + "," +
            self.chain_id + "," + self.i_code + ","
            for a in self.atoms:
                if a == None:
                    s += "N,"
                else:
                    s += a.to_str() + ","
        return s
