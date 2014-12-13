import numpy as np
from . import atom
from . import residue
from . import residue_type
from . import chain
from . import structure

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


def str_to_residue(s):
    """
    creates an residue from string generated from r.to_str()
    """
    spl = s.split(",")
    rtype = residue_type.get_rtype(spl[0])
    r = residue.Residue(rtype, spl[1], int(spl[2]), spl[3], spl[4])

    atoms = []
    for i in range(5, len(spl)-1):
        a = str_to_atom(spl[i])
        atoms.append(a)
    r.setup_atoms(atoms)
    return r


def str_to_chain(s):
    """
    creates a chain from string generated from chain.to_str()
    """
    spl = s.split(";")
    c = chain.Chain()
    residues = []
    for r_str in spl[:-1]:
        r = str_to_residue(r_str)
        residues.append(r)
    c.residues = residues
    return c


def str_to_structure(s):
    spl = s.split(":")
    struct = structure.Structure()
    chains = []
    for c_str in spl[:-1]:
        c = str_to_chain(c_str)
        chains.append(c)

    struct.chains = chains
    return struct

