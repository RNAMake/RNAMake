import numpy as np
from . import atom


def str_to_atom(s):
    """
    converts string to atom.Atom object format "AtomName X Y Z"
    :params s: string containing atom elements
    :type   s: string
    .. code-block:: python
        >>>str_to_atom("P 1.0 2.0 3.0")
        <Atom(name='P', coords='1.0 2.0 3.0')>
    """
    spl = s.split()
    coords = [float(x) for x in spl[1:]]
    return atom.Atom(spl[0], np.array(coords))
