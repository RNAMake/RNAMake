#coding=UTF-8

"""
This is the file to generate a tree storing the Fourier matrix element
of each member of an ensemble, or generating it if there is no enough RAM.
"""

from rnamake import base, option, motif_state_ensemble_tree,motif_state_tree
from rnamake import secondary_structure_parser, motif_type, motif_tree,exceptions
from rnamake import resource_manager as rm
from rnamake import se3util as se3

class MotifStateEnsembleConvolution:
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
    def __init__(self, mset, grid_size=(1,1,1,1,1,1),grid_unit=1):
        self.mset = mset
        self.mst_ = mset.to_mst()
        self.grid_size = grid_size
        self.grid_unit = grid_unit
        self.map = se3.SE3Map(self.grid_size, self.grid_unit)
        self.done_ = False
        self.ni_f_ = 1
        self.ni_l_ = mset.tree.last_node.index
        self.ei_f_ = 1
        self.ei_l_ = 1

    def run(self):

        # init
        self.map.place_ensemble(self.mset.get_node(self.ni_f_).data)

        #iter
        for en in self.mset.tree.nodes[2:]:#use iterator, start from one
            #if statement
            temp_map = se3.SE3Map(self.grid_size, self.grid_unit)
            temp_map.place_ensemble(en.data)
            # TODO order?
            self.map = temp_map * self.map
        self.done = True
    def get_energy(self):
        if not self.done:
            print("run me first")
            return 0
        self.map.get_state_probability(self.mst.get_node())
