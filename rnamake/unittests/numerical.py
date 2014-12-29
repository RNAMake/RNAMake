import rnamake.atom

def are_floats_equal(f1,f2):
    """
    test whether two floating point values are the same, cannot use ==,
    returns true if they are the same

    :param f1: first float
    :param f2: second float

    :type f1: float
    :type f2: float

	.. code-block:: python
        >>> are_floats_equal(1,1)
        1

        >>> are_floats_equal(1,1.00000000000001)
        1

        >>> are_floats_equal(1.001,1.001)
        1

        >>> are_floats_equal(1,1.00001)
        0
    """

    if abs(f1-f2) < 0.00001:
        return 1
    else:
        return 0


def are_points_equal(p1, p2):
    """
    test whether two points are equal, this is done by checking whether each
    element is the same between

    :param p1: point one
    :param p2: point two

    :type p1: list / np.array
    :type p2: list / np.array

	.. code-block:: python
        >>> are_points_equal([1,0,0],[1,0,0])
        1

        >>> are_points_equal([1,0,0],[1,0,1])
        0

    """
    for i in range(3):
        if not are_floats_equal(p1[i], p2[i]):
            return 0

    return 1


def are_matrices_equal(m1, m2):
    for i in range(len(m1)):
        for j in range(len(m1[i])):
            if not are_floats_equal(m1[i][j], m2[i][j]):
                return 0
    return 1


def are_atom_equal(a1, a2):
    return are_points_equal(a1.coords, a2.coords) and a1.name == a2.name
