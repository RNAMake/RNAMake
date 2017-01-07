import math
import numpy as np


def normalize(v):
    """
    gets a unit vector of the vector supplied

    :param v: the vector you wish to normalize
    :type v: np.array

    :return: a normalize vector
    :rtype: np.array
    """

    return v / math.sqrt( v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def distance(p1, p2):
    """
    returns the distance between two 3D points

    :param p1: point 1
    :param p2: point 2

    :type p1: np.array
    :type p2: np.array

    :return: distance betwene p1 and p2
    :rtype: float
    """
    return math.sqrt((p2[0] - p1[0]) ** 2 +
                     (p2[1] - p1[1]) ** 2 +
                     (p2[2] - p1[2]) ** 2)


def matrix_distance(r1, r2):
    """
    a quick way to evaluate how different two matrices are from each other.
    Calculates the sum of the absolute difference between each element
    between the matrices.

    :param r1: First 3x3 matrix
    :param r2: Second 3x3 matrix

    :type r1: 3x3 np.array
    :type r2: 3x3 np.array

    :return: difference between the matrices
    :rtype: float
    """

    dist = 0
    for i in range(3):
        for j in range(3):
            dist += abs(r1[i][j] - r2[i][j])
    return dist


def center(atoms):
    """
    returns the euclidean center of the coordinates of a list of atoms

    :params atoms: list of atoms to calculate center
    :type atoms: list of atoms

    :return: euclidean center of the atoms
    :rtype: np.array
    """

    center = np.array([0.0, 0.0, 0.0])
    count = 0
    for a in atoms:
        if a is None:
            continue
        center += a.coords
        count += 1

    return center / float(count)


def center_points(points):
    """
    returns the center of a list of numpy arrays

    :params points: list of points to calculate center
    :type points: list of np.arrays

    :return: euclidean center of the points
    :rtype: np.array
    """

    center = np.array([0.0, 0.0, 0.0])

    for p in points:
        center += p

    return center / float(len(points))


def base_dir(path):
    """
    gets the base directory that a directory or file is located

    :param path: full path of file or directory
    :type path: str

    :return: the base directory
    :rtype: str
    """

    path_spl = path.split("/")
    return "/".join(path_spl[:-1]) + "/"


def filename(path):
    """
    gets the filename of a path, excluding the rest of the path string

    :param path: the path you want the filename for
    :type path: str

    :return: the filename of a path
    :rtype: str
    """

    path_spl = path.split("/")
    return path_spl[-1]


def wc_bp(bp, s):
    """
    checks to see if the residues in a basepair match a watson-crick basepair

    :param bp: the basepair you wish to check to see if is watson-crick
    :type bp: basepair.Basepair

    :return: 1 if is watson-crick 0 if not
    :rtype: int
    """
    res1 = s.get_residue(uuid=bp.res1_uuid)
    res2 = s.get_residue(uuid=bp.res2_uuid)

    bp_str = res1.short_name() + res2.short_name()
    wc = "GC,CG,AU,UA".split(",")
    if bp_str in wc:
        return 1
    else:
        return 0


def gu_bp(bp, s):
    """
    checks to see if the residues in a basepair match a GU basepair

    :param bp: the basepair you wish to check to see if is GU
    :type bp: basepair.Basepair

    :return: 1 if is GU 0 if not
    :rtype: int
    """
    res1 = s.get_residue(uuid=bp.res1_uuid)
    res2 = s.get_residue(uuid=bp.res2_uuid)

    bp_str = res1.short_name() + res2.short_name()
    if bp_str == "GU" or bp_str == "UG":
        return 1
    else:
        return 0


def unitarize(R):
    """
    Enforce unitarity of the input matrix using Gram-Schmidt process.
    This function overwrites the input matrix.

    :param R: input matrix to be unitarized.
    :type  R: 3x3 np.array

    :return: Unitarized input matrix.
    :rtype:   3x3 np.array

    """
    if len(R.shape) == 2:  # 1D case
        R[0] /= math.sqrt(R[0].dot(R[0]))
        R[1] -= R[1].dot(R[0]) * R[0]
        R[1] /= math.sqrt(R[1].dot(R[1]))
        R[2] -= R[2].dot(R[0]) * R[0]
        R[2] -= R[2].dot(R[1]) * R[1]
        R[2] /= math.sqrt(R[2].dot(R[2]))
    else:  # 2D case
        R[:, 0] /= np.sqrt(np.sum(R[:, 0] ** 2, axis=1))[:, np.newaxis]
        R[:, 1] -= (
            np.einsum('ij,ij->i', R[:, 1], R[:, 0])[:, np.newaxis] * R[:, 0])
        R[:, 1] /= np.sqrt(np.sum(R[:, 1] ** 2, axis=1))[:, np.newaxis]
        R[:, 2] -= (
            np.einsum('ij,ij->i', R[:, 2], R[:, 0])[:, np.newaxis] * R[:, 0])
        R[:, 2] -= (
            np.einsum('ij,ij->i', R[:, 2], R[:, 1])[:, np.newaxis] * R[:, 1])
        R[:, 2] /= np.sqrt(np.sum(R[:, 2] ** 2, axis=1))[:, np.newaxis]
    return R



