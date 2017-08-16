#coding=UTF-8

"""
This is the file to generate a tree storing the Fourier matrix element
of each member of an ensemble, or generating it if there is no enough RAM.
"""

from rnamake import base, option, motif_state_ensemble_tree,motif_state_tree
from rnamake import secondary_structure_parser, motif_type, motif_tree,exceptions
from rnamake import resource_manager as rm
from rnamake import se3util as se3
import numpy as np

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
        self.ni_f_ = 2
        self.ni_l_ = mset.tree.last_node.index
        self.slidetrack = []

    def run(self):

        # init
        self.map.place_ensemble(self.mset.get_node(self.ni_f_).data)
        self.slidemap(self.ni_f_)

        #iter
        for en in self.mset.tree:#use iterator, start from one
            if en.index<=self.ni_f_:
                continue
            if en.index> self.ni_l_:
                self.done_ = True
                break
            #if statement
            temp_map = se3.SE3Map(self.grid_size, self.grid_unit)
            temp_map.place_ensemble(en.data)
            self.map = temp_map * self.map
            self.slidemap(en.index)
    def get_prob(self):
        if not self.done:
            print("run me first")
            return 0
        

    def slidemap(self,ni):
        gl = se3.MotifGaussianList(self.mset)
        meand = gl.mgl[ni].mean[:3,3]
        meangrid = np.around(meand/self.grid_unit)
        assert np.all(abs(meangrid)<self.grid_size)
        if meangrid[0]>0:
            self.map.data[:-meangrid[0],:,:,:,:,:] = self.map.data[meangrid[0]:,:,:,:,:,:]
            self.map.data[-meangrid[0]:,:,:,:,:,:] = 0
        if meangrid[0]<0:
            self.map.data[-meangrid[0]:,:,:,:,:,:] =  self.map.data[:meangrid[0],:,:,:,:,:]
            self.map.data[:-meangrid[0],:,:,:,:,:] = 0
        if meangrid[1]>0:
            self.map.data[:,:-meangrid[1],:,:,:,:] = self.map.data[:,meangrid[1]:,:,:,:,:]
            self.map.data[:,-meangrid[1]:,:,:,:,:] = 0
        if meangrid[1]<0:
            self.map.data[:,-meangrid[1]:,:,:,:,:] =  self.map.data[:,:meangrid[1],:,:,:,:]
            self.map.data[:,:-meangrid[1],:,:,:,:] = 0
        if meangrid[2]>0:
            self.map.data[:,:,:-meangrid[2],:,:,:] = self.map.data[:,:,meangrid[2]:,:,:,:]
            self.map.data[:,:,-meangrid[2]:,:,:,:] = 0
        if meangrid[2]<0:
            self.map.data[:,:,-meangrid[2]:,:,:,:] =  self.map.data[:,:,:meangrid[2],:,:,:]
            self.map.data[:,:,:-meangrid[2],:,:,:] = 0
        self.slidetrack.append(meangrid)



