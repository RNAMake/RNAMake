#coding=UTF-8

"""
This is the file to generate a tree storing the Fourier matrix element
of each member of an ensemble, or generating it if there is no enough RAM.
"""

from rnamake import base, option, motif_state_ensemble_tree
from rnamake import secondary_structure_parser, motif_type, motif_tree,exceptions
from rnamake import resource_manager as rm
import se3util as se3

def motif_state_ensemble_convolution(mset, grid_size=(1,1,1,1,1,1),grid_unit=1):
    """
    :param mset: MotifStateEnsembleTree to calculate
    :param grid_size: how many samples to sample along each axis
    :param grid_unit: what is the length between two grid points
    :type mset: motif_state_ensemble_tree.MotifStateEnsembleTree
    :type grid_size: iterable of int
    :type grid_unit: float
    :return: result distribution in ndarray
    :rtype: se3.SE3Map
    """

