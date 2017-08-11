import numpy as np

from rnamake import base, option, motif_state_ensemble_tree
from rnamake import motif_ensemble
from rnamake import motif,basepair
from rnamake import transformations
from rnamake import fconv3d

class SE3Map:
    """
    A map in SE(3) with all the information specified
    :param grid_size: How many grid along each axis
    :param grid_unit: How long is one grid
    :param data: numpy array actually storing the data
    :type grid_size: iterable of int
    :type grid_unit: float
    :type data: numpy.ndarray
    """
    def __init__(self, grid_size, grid_unit):
        self.grid_size = grid_size
        self.grid_unit = grid_unit
        self.data = np.zeros(grid_size)

    def __mul__(self, other):
        """
        operator reload for convolution
        :param other: right operand
        :type other: SE3Map
        :return:
        :rtype: SE3Map
        """
        assert self.grid_size == other.grid_size\
            and self.grid_unit == other.grid_unit
        res = SE3Map(self.grid_size, self.grid_unit)
        res.data = fconv3d.fconv3d_broadcast.c3d(self.data, other.data)
        return res


    #TODO
    def get_state_energy(self, basepair_state):
        """

        :param basepair_state:
        :type basepair_state: basepair.BasepairState
        :return:
        :rtype: float
        """
        pass




    def place_ensemble(self, mse):
        """
        :param mse: motif state ensemble to place in the map
        :type mse: motif_ensemble.MotifStateEnsemble
        :return: self
        """
        for mem in mse.members:
            self.place_motif_state(mem.motif_state,mem.energy)
        return self


    def place_motif_state(self, ms,nrg):
        """
        :param ms: motif state to place
        :param nrg: energy related to that MS
        :type ms: motif.MotifState
        :type nrg: float
        :return: None
        """
        # TODO align_to_org(ms)
        s = ms.end_states[0]
        grid_ndx = self.state_to_grid_ndx(s)

        probability = self.nrg_to_probability(nrg)

        self.data[grid_ndx] = probability


    #TODO
    def state_to_grid_ndx(self,basepair_state):
        """
        converting state to index on the grid
        :param basepair_state:
        :type basepair_state: basepair.BasepairState
        :return:
        :rtype: list of ints, shape(6,)
        """
        eu = transformations.euler_from_matrix(
            basepair_state.r,axes='szxz'
        )
        d = basepair_state.d
        a = self.grid_unit
        n = self.grid_size
        ndx = []
        return ndx


    #TODO
    def nrg_to_probability(self, energy):
        """
        conversion from energy to  probability density
        :param energy:
        :type energy: float
        :return: PD
        :rtype: float
        """
        pass