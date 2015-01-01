import math
import numpy as np

def normalize(v):
    return v / math.sqrt( v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

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


def matrix_distance(r1, r2):
    dist = 0
    for i in range(3):
        for j in range(3):
            dist += abs(r1[i][j] - r2[i][j])
    return dist

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


def base_dir(path):
    path_spl = path.split("/")
    return "/".join(path_spl[:-1]) + "/"


def filename(path):
    path_spl = path.split("/")
    return path_spl[-1]


def wc_bp(bp):
    bp_str = bp.res1.rtype.name[0] + bp.res2.rtype.name[0]
    wc = "GC,CG,AU,UA".split(",")
    if bp_str in wc:
        return 1
    else:
        return 0


def gu_bp(bp):
    bp_str = bp.res1.rtype.name[0] + bp.res2.rtype.name[0]
    if bp_str == "GU" or bp_str == "UC":
        return 1
    else:
        return 0

