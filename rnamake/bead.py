import util
import exceptions
import numpy as np

import basic_io

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

    @staticmethod
    def valid_bead_type(btype):
        if btype > 2 or btype < 0:
            return 0
        else:
            return 1


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

    __slots__ = ["__center", "__btype"]

    def __init__(self, center, btype):
        if not BeadType.valid_bead_type(btype):
            raise exceptions.ResidueException("not a valid BeadType value: " + str(btype))

        self.__center, self.__btype = center, btype

    @classmethod
    def from_str(cls, s):
        spl = s.split(",")
        center = basic_io.str_to_point(spl[0])
        return cls(center, int(spl[1]))

    @classmethod
    def copy(cls, b):
        """
        returns a deep copy of current bead object

        :returns: copy of bead object
        :rtype: Bead
        """

        return cls(np.copy(b.__center), b.__btype)

    def __repr__(self):
        center = basic_io.point_to_str(self.__center)
        return "<Bead(btype='%s', center='%s')>" % (self.type_name(), center)

    def type_name(self):
        """
        returns name of btype in string form

        :returns: name of btype type
        :rtype: str
        """
        if   self.__btype == 0:
            return "PHOSPHATE"
        elif self.__btype == 1:
            return "SUGAR"
        elif self.__btype == 2:
            return "BASE"
        else:
            raise exceptions.ResidueException("invalid bead type: " + str(self.__btype))

    def distance(self, b):
        return util.distance(self.__center, b.__center)

    def move(self, p):
        self.__center += p

    def transform(self, t):
        self.__center = np.dot(self.__center, t.rotation().T) + t.translation()

    def fast_transform(self, r, t):
        self.__center = np.dot(self.__center, r) + t

    def to_str(self):
        return basic_io.point_to_str(self.__center) + "," + str(self.__btype)

    @property
    def center(self):
        return np.copy(self.__center)

    @property
    def btype(self):
        return self.__btype
