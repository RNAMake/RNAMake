PDBLINE_GE100K = ('%-6s%5d %-4s%1s%-4s%1s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '%4s%2s\n')


def point_to_str(p):
    """
    :param p: point to print
    :type  p: list

    returns string of each element seperated by a space
    .. code-block:: python
        >>>p = [0,1,2]
        >>>point_to_str(p)
        0 1 2

    """
    return str(p[0]) + " " + str(p[1]) + " " + str(p[2])
