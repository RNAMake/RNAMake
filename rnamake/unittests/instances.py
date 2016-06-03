import numpy as np
from rnamake import settings, residue, chain, structure
from rnamake.unittests import files
import rnamake.structure
from rnamake import transform, transformations


def motif():
    path = rnamake.settings.MOTIF_DIRS + "base.motif"
    return rnamake.motif.file_to_motif(path, rts)


def residue():
    path = settings.UNITTEST_PATH + "resources/res_strs.dat"
    f = open(path)
    l = f.readline()
    f.close()

    return residue.Residue.from_str(l)


def chain():
    path = settings.UNITTEST_PATH + "resources/chain_strs.dat"
    f = open(path)
    l = f.readline()
    f.close()

    return chain.Chain.from_str(l)


def structure():
    s = rnamake.structure.structure_from_pdb(files.P4P6_PDB_PATH)
    return s


def transform_indentity():
    t = transform.Transform()
    return t


def transform_random():
    r = transformations.random_rotation_matrix()[:3, :3]
    d = np.random.random([3])*50
    return transform.Transform(r, d)


