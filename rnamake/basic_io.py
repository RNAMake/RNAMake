import numpy as np

PDBLINE_GE100K = ('%-6s%5d %-4s%1s%-4s%1s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '%4s%2s\n')


def point_to_str(p):
    return str(p[0]) + " " + str(p[1]) + " " + str(p[2])
