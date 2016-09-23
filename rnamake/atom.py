import basic_io
import numpy as np


class Atom(object):
    """Stores atomic information from pdb file, design is to be extremely
    lightweight only storing the atom name and coordinates.

    :param name: name of atom
    :param coords: 3d coordinates of atom's position

    :type name: str
    :type coords: numpy.array


    :attributes:

    `name` : str
        Atomic name
    `coords` : np.array
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
    __slots__ = ["name", "coords"]

    def __init__(self, name, coords):
        """returns new atom.Atom object"""

        self.name, self.coords = name, coords

    def __repr__(self):
        """returns string representation of object"""

        coords = basic_io.point_to_str(self.coords)
        return "<Atom(name='%s', coords='%s')>" % (self.name, coords)

    def copy(self):
        """Deep copies the current atom instance.

        :returns: an Atom object

        :examples:

        .. code-block:: python

            >>> a = Atom("P",[1.0,2.0,3.0])
            >>> a_copy = a.copy()
            >>> print a_copy.name
            P

        """

        coords = np.array(self.coords)
        return Atom(self.name, coords)

    def to_str(self):
        """returns string version of atom.

        :returns: str

        :examples:

        .. code-block:: python

            >>> atom = atom.Atom("H1", np.array([0, 1, 2]))
            >>> string = atom.to_str()
            "H1 0.0 1.0 2.0"
        """
        return self.name + " " + basic_io.point_to_str(self.coords)

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
                self.coords[0],
                self.coords[1],
                self.coords[2])
        return s
