import numpy as np

PDBLINE_GE100K = ('%-6s%5d %-4s%1s%-4s%1s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '%4s%2s\n')


def point_to_str(p):
    """
    converts numpy array or list into string

    :param p: point to print
    :type  p: list

    returns string of each element seperated by a space
    .. code-block:: python
        >>>p = [0,1,2]
        >>>point_to_str(p)
        0 1 2

    """
    return str(p[0]) + " " + str(p[1]) + " " + str(p[2])


def points_to_str(points):
    """
    converts a list of points into a string seperated by spaces

    .. code-block:: python
        >>> points = [[1, 0, 1], [0, 0, 0]]
        >>> points_to_str(points)
        1 0 1 0 0 0

    """
    s = ""
    for p in points:
        s += point_to_str(p) + " "
    return s


def matrix_to_str(m):
    """
    converts numpy array that is multidimensional to a string, assumes that
    each column is length 3, this might be an issue later

    :params

    .. code-block:: python
        >>> m = np.eye(3)
        >>> matrix_to_str(m)
        1 0 0 0 1 0 0 0 1

    """
    s = ""
    for i in range(len(m)):
        for j in range(3):
            s += str(m[i][j]) + " "
    return s


def str_to_point(s):
    p = np.array([float(x) for x in s.split()])
    return p


def str_to_matrix(s):
    m = np.array([float(x) for x in s.split()])
    rows = int(len(m)/3)
    m = m.reshape(rows,3)
    return m

