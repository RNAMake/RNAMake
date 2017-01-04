import basic_io
import exceptions

import numpy as np


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

            >>> str_to_atom("P 1.0 2.0 3.0")
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


    def to_str(self):
        """returns string version of atom.

        :returns: str

        :examples:

        .. code-block:: python

            >>> atom = atom.Atom("H1", np.array([0, 1, 2]))
            >>> string = atom.to_str()
            "H1 0.0 1.0 2.0"
        """
        return self.__name + " " + basic_io.point_to_str(self.__coords)

    def to_pdb_str(self, acount=1):
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

    def transform(self, t):
        self.__coords = np.dot(self.__coords, t.rotation().T) + t.translation()

    def move(self, p):
        self.__coords += p

    @property
    def coords(self):
        return self.__coords

    @property
    def name(self):
        return self.__name

    @coords.setter
    def coords(self, c):
        raise exceptions.AtomException(
            "cannot set coords externally!, either use move or transform, "
            "this is to keep encapsulation")

    @name.setter
    def name(self, n):
        raise exceptions.AtomException(
            "cannot set name externally!, this is to keep encapsulation")

