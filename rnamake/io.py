import numpy as np
from . import atom

def str_to_atom(s):
    """
    """

    spl = s.split()
    coords = [float(x) for x in spl[1:]]
    return atom.Atom(spl[0], np.array(coords))

