import math

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
