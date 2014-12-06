import math
import numpy as np


def distance(p1, p2):
    """
    returns the distance between two 3D points

    :param p1: point 1
    :param p2: point 2
    :type p1: list/array
    :type p2: list/array
    """
    return math.sqrt((p2[0] - p1[0]) ** 2 +
                     (p2[1] - p1[1]) ** 2 +
                     (p2[2] - p1[2]) ** 2)


def center(atoms):
    """
    returns the center of a list of atoms

    :params atoms: list of atoms to calculate center
    :type atoms: list of atoms
    """

    center = np.array([0.0, 0.0, 0.0])

    for a in atoms:
        if a is None:
            continue
        center += a.coords

    return center / float(len(atoms))
