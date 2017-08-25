import numpy as np

from rnamake import base, option, motif_state_ensemble_tree
from rnamake import motif_ensemble
from rnamake import motif,basepair,transform
from rnamake import transformations
from rnamake import fconv3d
from rnamake.fconv3d import auxi
from rnamake.fconv3d import envis as ev
from numpy import linalg as la
import copy
from numba import jit

class SE3Map(object):
    """
    A map in SE(3) with all the information specified
    :param grid_size: How many grid along each axis
    :param grid_unit: How long is one grid
    :param data: numpy array actually storing the data
    :type grid_sizes: iterable of int,length 4
    :type grid_unit: float
    :type data: numpy.ndarray
    """
    def __init__(self, grid_sizes, grid_unit):
        self.grid_sizes = grid_sizes
        self.grid_unit = grid_unit
        self.data = np.zeros([grid_sizes[0], grid_sizes[0], grid_sizes[0],
                              grid_sizes[1], grid_sizes[2], grid_sizes[3]])
        org = (self.grid_sizes[0] - 1) / 2
        auxi.org = [org,org,org]
        self.gen_be_weight()

    def __mul__(self, other):
        """
        operator reload for convolution
        :param other: right operand
        :type other: SE3Map
        :return:
        :rtype: SE3Map
        """
        assert self.grid_sizes == other.grid_sizes\
            and self.grid_unit == other.grid_unit
        res = SE3Map(self.grid_sizes, self.grid_unit)
        res.data = fconv3d.fconv3d.c3d(self.data, other.data)
        res.data /= self.nmlz()
        res.data[res.data<0]=0
        res.data[np.isclose(res.data,0)] = 0
        assert not np.all(res.data == 0) # There should still be some data
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
        netc=0
        for mem in mse.members:
            netc += mem.count
        for mem in mse.members:
            probability = float(mem.count)/netc
            self.place_motif_state(mem.motif_state,probability)
        ev.envis6d(self.data)
        return self

    #finished
    def place_motif_state(self, ms, prob):
        """
        :param ms: motif state to place
        :param prob: probability related to that MS
        :type ms: motif.MotifState
        :type prob: float
        :return: None
        """
        start_matrix = state_to_matrix(ms.end_states[0])
        final_matrix  = state_to_matrix( ms.end_states[1])
        step_matrix = la.inv(start_matrix).dot(final_matrix)
        # probability = nrg_to_probability(prob)
        # e_r, e_d = ms.end_states[1].r, ms.end_states[1].d
        # s_r, s_d = ms.end_states[0].r, ms.end_states[0].d
        # t = transform.Transform(s_r.T, )
        # s = ms.end_states[0]
        grid_ndx = self.matrix_to_grid_ndx(step_matrix)
        # probability = nrg_to_probability(prob)
        self.data[grid_ndx] += prob/self.voxel(grid_ndx)
        # print 'motif state to place\n',step_matrix,'\n'
        # print 'it is placed at', grid_ndx



    #finished
    def matrix_to_grid_ndx(self, state_matrix):
        """
        converting state to index on the grid
        :param state_matrix:
        :type state_matrix: numpy.ndarray
        :return:
        :rtype: tuple of ints, shape(6,)
        """
        a = self.grid_unit
        n = np.asarray(self.grid_sizes).astype('float')
        eu = transformations.euler_from_matrix(
            state_matrix[:3,:3],axes='szxz'
        )[::-1]
        d = state_matrix[:3,3].astype('float')
        nd = np.around(d/a)
        if not np.all(abs(nd)<=(n[0]-1)/2):
            print 'shift index required',nd
            assert 0 # keep index in range
        nd += (n[0]-1)/2

        neu = np.zeros(3)
        neu[0]=np.around(eu[0]*n[1]/2/np.pi)
        if neu[0]==n[1]:
            neu[0] = 0
        neu[1]=int(eu[1]*n[2]/np.pi)
        if neu[1]==n[2]:
            print 'beta out of index!'
            neu[1] = n[2] -1
        neu[2] = np.around(eu[2]*n[3]/2/np.pi)
        if neu[2] == n[3]:
            neu[2] = 0

        return (nd[0],nd[1],nd[2],neu[0],neu[1],neu[2])
    #finished
    def voxel(self,grid_ndx):
        """
        :type grid_ndx: tuple of ints, shape(6,)
        :param grid_ndx:
        :return:
        """
        beta_grid_size = self.grid_sizes[2]
        if grid_ndx[4] == 0:
            beta_w = 1-np.cos(np.pi/2/beta_grid_size)
        elif grid_ndx[4] == beta_grid_size-1:
            beta_w = np.cos(np.pi*(beta_grid_size-3.0/2)/beta_grid_size)+1
        else:
            beta_w = np.cos(np.pi/beta_grid_size*(grid_ndx[4]-0.5))-np.cos(np.pi/beta_grid_size*(grid_ndx[4]+0.5))
        res = beta_w * self.grid_unit**3
        return res

    def gen_be_weight(self):
        beta_grid_size = self.grid_sizes[2]
        self.bwt = np.zeros(beta_grid_size,dtype='float32')
        for be in range(beta_grid_size):
            if be == 0:
                self.bwt[be] = 1 - np.cos(np.pi / 2 / beta_grid_size)
            elif be == beta_grid_size - 1:
                self.bwt[be] = np.cos(np.pi * (beta_grid_size - 3.0 / 2) / beta_grid_size) + 1
            else:
                self.bwt[be] = np.cos(np.pi / beta_grid_size * (be - 0.5)) - np.cos(
                    np.pi / beta_grid_size * (be + 0.5))



    def trans_reg(self):
        # regression to translational distributions
        return np.sum(self.data*self.bwt[np.newaxis,np.newaxis,np.newaxis,np.newaxis,:,np.newaxis],
                      axis=(3,4,5))



    def nmlz(self):
        res = 0
        # for x1 in range(self.grid_sizes[0]):
        #     for x2 in range(self.grid_sizes[0]):
        #         for x3 in range(self.grid_sizes[0]):
        for al in range(self.grid_sizes[1]):
            for be in range(self.grid_sizes[2]):
                for ga in range(self.grid_sizes[3]):
                    res += np.sum(self.data[:,:,:,al,be,ga]*self.voxel([0,0,0,al,be,ga]))
        return res
        pass


class MotifGaussianList(object): # finished
    __slots__ = ['mset', 'mgl']



    def __init__(self,mset):
        self.set_from_mset(mset)


    def set_from_mset(self,mset):
        """

        :param mset:
        :type mset: motif_state_ensemble_tree.MotifStateEnsembleTree
        :return:
        """
        self.mgl=[]
        for en in mset:
            states = np.zeros([4,4,len(en.data.members)])
            counts = np.zeros(len(en.data.members))
            for i,msm in enumerate(en.data.members):
                final_state= state_to_matrix(msm.motif_state.end_states[1])
                start_state= state_to_matrix(msm.motif_state.end_states[0])
                states[:,:,i] = np.dot(la.inv(start_state),final_state)
                counts[i] = msm.count
            # print en.index,'\t',np.sum(counts)
            self.mgl.append(self.mg_from_sc(states,counts))


    def mg_from_sc(self, st, ct):
        # print ct
        assert type(st) == np.ndarray \
               and st.shape[:2] ==(4,4)  \
               and ct.ndim == 1
        mean = np.average(st, axis=2,weights=ct)
        # mean = st[:,:,0]
        mean[:3,:3] /= (la.det(mean[:3,:3])**(1.0/3))
        covar = np.cov(matrix_to_chi(st-mean[:,:,np.newaxis]), fweights=ct)
        if not np.all(np.isfinite(covar)):
            assert np.all(np.isnan(covar))
        return MotifGaussian(covar,mean)


    def get_mg(self,ni1, ni2):
        res_mg = copy.copy(self.mgl[ni1])
        for i in range(ni1+1,ni2+1):
            res_mg = res_mg*self.mgl[i]
            # print '|SIGMA| at the end of step %d = '%i,la.det(res_mg.SIGMA)
        return res_mg



class MotifGaussian(object): # finished
    __slots__ = ['SIGMA','mean']

    def __init__(self,SIGMA,mean):
        assert (type(SIGMA), type(mean) == np.ndarray, np.ndarray)\
            and (SIGMA.shape, mean.shape == (6,6),(4,4))
        self.SIGMA=np.nan_to_num(SIGMA).copy()
        self.mean = mean.copy()

    def __copy__(self):
        print 'MotifGaussian copied'
        return MotifGaussian(self.SIGMA,self.mean)

    def __mul__(self,other):
        """
        :type other: MotifGaussian
        :param other:
        :return:
        """
        assert self.mean.shape ==(4,4)
        res_mean = np.dot(self.mean,other.mean)
        g2_inv = la.inv(other.mean)
        r = g2_inv[:3,:3]
        d = g2_inv[:3, 3]
        ad = np.zeros([6,6])
        ad[:3,:3] = r
        ad[3:,3:] = r
        ad[3:,:3] = np.cross(d,r,axisb=0)
        res_cov = np.dot(np.dot(ad,self.SIGMA),ad.T)+other.SIGMA
        return MotifGaussian(res_cov,res_mean)


    def eval(self,chi):
        assert type(chi) == np.ndarray\
            and chi.shape == (6,)
        chi = chi.astype('float')

        return 1.0/((2*np.pi)**3*np.sqrt(la.det(self.SIGMA)))*\
               np.exp(-1.0/2*np.dot(np.dot(chi.T,la.inv(self.SIGMA)),chi))#[0,0]

    def eval_double(self,chi,TAU):
        assert type(chi) == np.ndarray\
        and chi.shape == (6,)\
        and type(TAU) == np.ndarray\
        and TAU.shape == (6,6)
        TAU = TAU.astype('float')
        A = la.inv(self.SIGMA)
        return 1.0/np.sqrt(la.det(TAU+A)*la.det(self.SIGMA))*\
                np.exp(-1.0/2*np.dot(np.dot(chi,TAU),chi))*np.exp(1.0/2*(np.dot(np.dot(np.dot(TAU,chi),
                                                                                       la.inv(A+TAU)),np.dot(TAU,chi))))




def nrg_to_probability(energy):
    """
    conversion from energy to  probability density
    :param energy:
    :type energy: float
    :return: PD
    :rtype: float
    """

    '''
    The following constant is copied from rnamake/setup/build_sqlite_libraries:337.
    Please check consistency if necessary.
    '''
    kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
    kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

    return np.exp(energy/(-kBT))

def prob_to_nrg(prob):
    kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
    kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)
    return kBT*np.log(prob)*6.02/100



def matrix_to_chi(t):
    chi1 = (t[2, 1] - t[1, 2]) / 2
    chi2 = (t[0, 2] - t[2, 0]) / 2
    chi3 = (t[1, 0] - t[0, 1]) / 2
    chi4, chi5, chi6 = t[0,3], t[1,3], t[2,3]
    chi =np.array([chi1,chi2,chi3,chi4,chi5,chi6])
    return chi

def state_to_matrix(bs):
    """
    :type bs: basepair.BasepairState
    :param bs:
    :return:
    """
    res = np.zeros([4,4]).astype('float')
    res[:3,:3] = bs.r.copy().T
    res[:3,3] = bs.d.copy()
    res[:3, :3] /= (la.det(res[:3, :3]) ** (1.0 / 3))
    res[3,3] = 1
    return res

def test_mul():
    from rnamake import unittests
    from rnamake.unittests import instances as inst
    motif1 = inst.motif()
    motif2 = inst.motif()
    st1 = motif1.ends[1].state()
    st2 = motif2.ends[1].state()
    ost2 = copy.copy(motif2.ends[0].state())
    st2c = copy.copy(st2)
    st1m = state_to_matrix(st1)
    st2m = state_to_matrix(st2)
    stm = np.dot(st1m,st2m)
    stminv = np.dot(st2m,st1m)
    motif.align_motif(motif1.ends[1].state(),motif2.ends[0],motif2)
    st5 = motif1.ends[0].state()
    st4 = motif1.ends[1].state()
    st3 = motif2.ends[0].state()
    st6 = motif2.ends[1].state()
    print 'st5',st5
    print 'st4',st4
    print 'st3',st3
    print 'st6',st6
    print 'st2c',st2c
    print 'ost2',ost2
    print 'st1m\n',st1m,'\n'
    print 'stm \n',stm, '\n'
    print 'stminv\n',stminv,'\n'
    motif1.to_pdb('test_motif1.pdb')
    motif2.to_pdb('test_motif2.pdb')

if __name__ == "__main__":
    test_mul()