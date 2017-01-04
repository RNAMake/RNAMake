import numpy as np
from rnamake import transform, transformations


def transform_indentity():
    t = transform.Transform()
    return t


def transform_random():
    r = transformations.random_rotation_matrix()[:3, :3]
    d = np.random.random([3])*50
    return transform.Transform(r, d)
