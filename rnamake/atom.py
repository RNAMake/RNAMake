from . import basic_io


class Atom(object):

    """
    stores atomic information from pdb file, design is to be extremely
    lightweight only storing the atom name and coordinates.

    :param name: atomic name
    :type name: str
    :param coords: atomic coordinates
    :type coords: list

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
        self.name, self.coords = name, coords

    def to_str(self):
        return self.name + " " + basic_io.point_to_str(self.coords)
