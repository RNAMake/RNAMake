import basic_io
import numpy as np

class Atom(object):

    """
    stores atomic information from pdb file, design is to be extremely
    lightweight only storing the atom name and coordinates.

    :param name: atomic name
    :param coords: atomic coordinates

    :type name : str
    :type coords : list

    Attributes
    ----------
    `name` : str
        Atomic name
    `coords` : Numpy array
        Atomic coordinates

    Examples

    .. code-block:: python
        >>>a = Atom("P",[1.0,2.0,3.0])
        >>>a.name
        P

        >>>a.coords
        [1.0 2.0 3.0]

        >>>print a
        <Atom(name ='P', coords='1.0 2.0 3.0')>

    """
    __slots__ = ["name", "coords"]

    def __init__(self, name, coords):
        """
        returns new rnamake.atom.Atom object
        """
        self.name, self.coords = name, coords

    def __repr__(self):
        coords = basic_io.point_to_str(self.coords)
        return "<Atom(name='%s', coords='%s')>" % (self.name, coords)

    def copy(self):
        coords = np.array(self.coords)
        return Atom(self.name, coords)

    def to_str(self):
        """
        returns string version of atom
        .. code-block:: python
            >>>atom = rnamake.atom.Atom("H1", np.array([0, 1, 2]))
            >>>string = atom.to_str()
            "H1 0.0 1.0 2.0"
        """
        return self.name + " " + basic_io.point_to_str(self.coords)

    def to_pdb_str(self, acount=1):
        """
        prints the current atom to a pdb string
        :params acount: the atom number of the atom in pdb file
        :type   acount: int
        .. code-block:: python
            >>>a = Atom("P",[1.0,2.0,3.0])
            >>>atom.to_pdb_str()
            ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P

            >>>atom.to_pdb_str(10)
            ATOM     10  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P
        """
        if self is None:
            return ""
        string = "ATOM {:6d}  P   C   A   1 {:11.3f}{:8.3f}{:8.3f}  1.00 62.18           P\n".format(
            acount,
            self.coords[0],
            self.coords[1],
            self.coords[2])
        return string
