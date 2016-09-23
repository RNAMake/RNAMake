import numpy as np
from rnamake import settings, io
from rnamake.unittests import files
import rnamake.structure
import rnamake.motif
from rnamake import transform, transformations


def motif():
    path = rnamake.settings.MOTIF_DIRS + "base.motif"
    return rnamake.motif.file_to_motif(path)


def residue():
    path = settings.UNITTEST_PATH + "resources/res_strs.dat"
    f = open(path)
    l = f.readline()
    f.close()

    return io.str_to_residue(l)


def chain():
    path = settings.UNITTEST_PATH + "resources/chain_strs.dat"
    f = open(path)
    l = f.readline()
    f.close()

    return io.str_to_chain(l)


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


def basepair():
    return motif().ends[0]


def baspairstate():
    return basepair().state()


def basepairstate_random():
    bp_state = baspairstate()
    t = transform_random()
    new_r, new_d, new_sug = bp_state.get_transformed_state(t.rotation(), t.translation())
    bp_state.set(new_r, new_d, new_sug)

    return bp_state


def secondary_structure_motif():
    return motif().secondary_structure


def secondary_structure_chain():
    return secondary_structure_motif().chains()[0]
