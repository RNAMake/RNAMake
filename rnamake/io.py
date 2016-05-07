import numpy as np

import atom
import residue
import residue_type
import chain
import structure

# TODO going to phase out this module in general and move each function to its
# respective module with the class its creating

def str_to_atom(s):
    """
    converts string to atom.Atom object format "AtomName X Y Z"

    :params s: string containing atom elements
    :type s: str

    .. code-block:: python

        >>> str_to_atom("P 1.0 2.0 3.0")
        <Atom(name='P', coords='1.0 2.0 3.0')>
    """
    spl = s.split()
    coords = [float(x) for x in spl[1:]]
    return atom.Atom(spl[0], np.array(coords))


def str_to_residue(s):
    """
    creates an residue from string generated from
    :func:`rnamake.residue.Residue.to_str`

    :param s: string containing stringifed residue
    :type s: str

    :returns: unstringifed residue object
    :rtype: residue.Residue

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
    creates an chain from string generated from
    :func:`rnamake.chain.Chain.to_str`

    :param s: string containing stringifed chain
    :type s: str

    :returns: unstringifed chain object
    :rtype: chain.Chain
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
    """
    creates an structure from string generated from
    :func:`rnamake.structure.Structure.to_str`

    :param s: string containing stringifed structure
    :type s: str

    :returns: unstringifed structure object
    :rtype: structure.Structure
    """

    spl = s.split(":")
    struct = structure.Structure()
    chains = []
    for c_str in spl[:-1]:
        c = str_to_chain(c_str)
        chains.append(c)

    struct.chains = chains
    return struct

