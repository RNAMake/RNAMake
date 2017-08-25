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
import copy
from scipy import linalg as la
from rnamake import transformations as tf
# from numba import jit

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
    def __init__(self, mset,ni1,ni2, grid_size=(1,1,1,1),grid_unit=1):
        self.mset = mset
        self.mst_ = mset.to_mst()
        self.grid_size = grid_size
        self.grid_unit = grid_unit
        self.map = se3.SE3Map(self.grid_size, self.grid_unit)
        self.done_ = False
        self.ni_f_ = ni1
        self.ni_l_ = ni2
        self.slidetrack = []
        self.gl = se3.MotifGaussianList(self.mset)
        self.gingers=[12] #tetraloops with only one member can be treated specially

    def run(self):


        # init
        self.map.place_ensemble(self.mset.get_node(self.ni_f_).data)
        # self.seemap()
        self.slidemap(self.ni_f_)
        # self.seemap()

        #iter
        for en in self.mset.tree.nodes:#use iterator, start from one
            print 'currently processing node # %d'%en.index
            if en.index<=self.ni_f_:
                continue
            if en.index in self.gingers:
                print 'tackling ginger %d naming %s'%(en.index,en.data.members[0].motif_state.name)
                self.dashmap(en.index)# may or may not be erratic
                continue
                #if statement
            temp_map = se3.SE3Map(self.grid_size, self.grid_unit)
            temp_map.place_ensemble(en.data)
            self.map = self.map*temp_map
            if not en.index == self.ni_l_:
                self.slidemap(en.index)
            # del temp_map
            assert np.all(np.isfinite(self.map.data))
            if en.index>= self.ni_l_:
                self.done_ = True
                lastdata = self.map.data
                slidetotal = np.zeros(3)
                for sld in self.slidetrack:
                    slidetotal += sld
                np.savez('run4',lastdata=lastdata,slidetotal=slidetotal)
                print 'self.done should be updated!!!'
                break
    def get_prob(self,mtx):
        m = mtx.copy()
        if not self.done_:
            print("run me first")
            return 0
        slidetotal = np.zeros(3)
        for sld in self.slidetrack:
            slidetotal += sld
        m[:3,3] -= slidetotal*self.grid_unit
        print 'slidetotal',slidetotal
        print m[:3,3]
        x1,x2,x3,al,be,ga = self.map.matrix_to_grid_ndx(m)
        netp =0
        # for x1p in np.arange(max(x1-2,0),min(x1+2+1,self.grid_size[0])):
            # for x2p in np.arange(max(x2-2,0),min(x2+2+1,self.grid_size[0])):
            #     for x3p in np.arange(max(x3 - 2, 0), min(x3 + 2 + 1, self.grid_size[0])):
                #     for alp in np.arange(al-2,al+2+1):
                #         for bep in np.arange(max(be-2,0),min(be+2,self.grid_size[2])):
                #             for gap in np.arange(ga-2,ga+2+1):
                #                 netp += self.map.data[x1p,x2p,x3p,alp,bep,gap]*self.map.voxel([x1p,x2p,x3p,alp,bep,gap])
        netp += self.map.data[x1, x2, x3, al, be, ga] * self.map.voxel([x1, x2, x3, al, be, ga])
        return netp/self.map.nmlz()


    def get_trans_pdf(self,mtx):
        m = mtx.copy()
        if not self.done_:
            print("run me first")
            return 0
        slidetotal = np.zeros(3)
        for sld in self.slidetrack:
            slidetotal += sld
        m[:3, 3] -= slidetotal * self.grid_unit
        print 'slidetotal', slidetotal
        print m[:3, 3]
        return self.map.trans_reg()[self.map.matrix_to_grid_ndx(m)]/self.map.nmlz()


    def slidemap(self,index):
        ni = copy.copy(index)
        meand_relative = self.gl.mgl[ni].mean[:3,3]
        meanr_ni = self.gl.get_mg(self.ni_f_,ni-1).mean[:3,:3]
        meandrot = np.dot(meanr_ni,meand_relative)
        if self.ni_l_-index< 8:
            end_state_mst = se3.state_to_matrix(self.mst_.get_node(1).data.cur_state.end_states[1])
            start_state_mst = se3.state_to_matrix(self.mst_.get_node(2).data.cur_state.end_states[0])
            step_mst = np.dot(la.inv(start_state_mst),end_state_mst)
            correction = step_mst[:3,3]-self.gl.get_mg(self.ni_f_,self.ni_l_).mean[:3,3]
            meandrot += correction/8.0
        # meand = self.mst_.get_node(ni).data.cur_state.end_states[1].d-\
        #     self.mst_.get_node(ni).data.cur_state.end_states[0].d
        # meand = self.mst_.get_node(ni).data.ref_state.end_states[1].d
        meangrid = np.around(meandrot/self.grid_unit)
        print 'mean ni1 to %d'%ni,self.gl.get_mg(self.ni_f_,ni).mean
        # print 'meand',meand
        print 'meandrot',meandrot
        # print 'meand_relative',meand_relative
        assert np.all(abs(meangrid)<self.grid_size[0])
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
        print 'meangrid',meangrid

    # @jit
    def dashmap(self,index):
        rang_x= self.map.data.shape[0]
        rang_al, rang_be, rang_ga = self.map.data.shape[3:]
        newmap = np.zeros(self.map.data.shape,dtype='complex64')
        dash_matrix = self.gl.mgl[index].mean
        mean_start = self.gl.get_mg(self.ni_f_,index-1).mean
        mean_end = self.gl.get_mg(self.ni_f_,index).mean
        mean_step = np.dot(la.inv(mean_start),mean_end)
        # print 'mean_step',mean_step
        # print 'dash_mat',dash_matrix
        # print 'mean_diff_d',mean_end[:3,3]-mean_start[:3,3]
        meand_relative = self.gl.mgl[index].mean[:3, 3]
        meanr_ni = self.gl.get_mg(self.ni_f_, index - 1).mean[:3, :3]
        meand = np.dot(meanr_ni, meand_relative)
        print 'meand',meand
        meangrid = np.around(meand / self.grid_unit)
        dash_slide = meangrid
        # print 'dash_slide',dash_slide
        # print 'mean_index-1.dot.dash_matrix diff d',mean_start.dot(dash_matrix)-mean_start
        self.slidetrack.append(meangrid)
        print 'meangrid',meangrid
        # print 'data maxima before dashing',self.map.data.max()
        loss = 0
        gain = 0
        for x1 in range(rang_x):
            for x2 in range(rang_x):
                for x3 in range(rang_x):
                    for al in range(rang_al):
                        for be in range(rang_be):
                            for ga in range(rang_ga):
                                data = self.map.data[x1,x2,x3,al,be,ga]
                                if np.isclose(data,0):
                                    continue
                                ang_al = al*2*np.pi/rang_al
                                ang_be = (be+0.5)*np.pi/rang_be
                                ang_ga = ga*2*np.pi/rang_ga
                                point_r = tf.euler_matrix(ang_ga,ang_be,ang_al,'szxz')[:3,:3]
                                point_d = np.array([x1,x2,x3])*self.grid_unit
                                point_mat = np.zeros([4,4],dtype='float32')
                                point_mat[:3,:3] = point_r
                                point_mat[:3,3] =point_d
                                new_mat = np.dot(point_mat,dash_matrix)
                                new_r = new_mat[:3,:3]
                                new_d = new_mat[:3,3]
                                # print 'actual diff d',new_d-point_d
                                new_eu = tf.euler_from_matrix(new_r,'szxz')[::-1]
                                new_eu_grid = np.zeros(3,dtype='float32')
                                new_eu_grid[0] = round(new_eu[0] * rang_al/2/np.pi)%rang_al
                                new_eu_grid[1] = int(new_eu[1]*rang_be/np.pi)
                                new_eu_grid[2] = round(new_eu[2]*rang_ga/2/np.pi)%rang_ga
                                new_d_grid = np.around(new_d/self.grid_unit)-dash_slide
                                if np.any(new_d_grid<0) or np.any(new_d_grid>=rang_x):
                                    loss+=1
                                    print 'out of range! the new grid'
                                    # print 'new_d',new_d
                                    # print 'new_grid',new_d_grid
                                    # print 'meangrid',meangrid
                                    # print 'meand relative',meand_relative
                                    # print 'mean d', self.gl.get_mg(self.ni_f_,index).mean
                                    # print 'point mat',point_mat
                                    continue
                                newmap[new_d_grid[0],new_d_grid[1],new_d_grid[2],
                                new_eu_grid[0],new_eu_grid[1],new_eu_grid[2]] = data
                                gain += 1
                                print '***normal data point in loop***'
        assert np.all(np.isfinite(newmap))
        # # assert np.sum(newmap)>0
        self.map.data = newmap
        print '\nloss %d versus gain %d\n'%(loss,gain)

    def seemap(self):
        print '--------SEEMAP--------'
        rang_x = self.map.data.shape[0]
        rang_al, rang_be, rang_ga = self.map.data.shape[3:]
        for x1 in range(rang_x):
            for x2 in range(rang_x):
                for x3 in range(rang_x):
                    for al in range(rang_al):
                        for be in range(rang_be):
                            for ga in range(rang_ga):
                                if np.isclose(self.map.data[x1,x2,x3,al,be,ga],0):
                                    continue
                                print 'seemap,',[x1,x2,x3,al,be,ga]

        print'--------END--------'

    def seetarget(self):
        from matplotlib import pyplot as plt

    def get_prob_analysis(self, mtx):
        m = mtx.copy()
        slidetotal = np.zeros(3)
        for sld in self.slidetrack:
            slidetotal += sld
        m[:3, 3] -= slidetotal * self.grid_unit
        print 'slidetotal', slidetotal
        print m[:3, 3]
        x1, x2, x3, al, be, ga = self.map.matrix_to_grid_ndx(m)
        print 'requested array index:', x1, x2, x3, al, be, ga
        netp = 0
        # for x1p in np.arange(max(x1-3,0),min(x1+3+1,self.grid_size[0])):
        #     for x2p in np.arange(max(x2-3,0),min(x2+3+1,self.grid_size[0])):
        #         for x3p in np.arange(max(x3 - 3, 0), min(x3 + 3 + 1, self.grid_size[0])):
        for alp in np.arange(self.grid_size[1]):
            for bep in np.arange(max(be-2,0),min(be+2,self.grid_size[2])):
                for gap in np.arange(al+ga-alp-1,al+ga-alp+1+1):
                    netp += self.map.data[x1,x2,x3,alp,bep,gap]*self.map.voxel([x1,x2,x3,alp,bep,gap])
                    # netp += self.map.data[x1p,x2p,x3p,alp,bep,gap]*self.map.voxel([x1p,x2p,x3p,alp,bep,gap])
        # netp += self.map.data[x1, x2, x3, al, be, ga] * self.map.voxel([x1, x2, x3, al, be, ga])
                    # netp += self.map.data[x1p, x2p, x3p, al, be, ga] * self.map.voxel([x1p, x2p, x3p, al, be, ga])
        print 'normalization factor:',self.map.nmlz()
        return netp





from rnamake import base, option, motif_state_ensemble_tree
from rnamake import secondary_structure_parser, motif_type, motif_tree, exceptions
from rnamake import resource_manager as rm

import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-fseq')
    parser.add_argument('-fss')
    parser.add_argument('-cseq')
    parser.add_argument('-css')

    args = parser.parse_args()
    return args

class SimulateTectos(base.Base):
    def __init__(self, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        #self.mset = self._get_mset()

    def setup_options_and_constraints(self):
        options = { 'fseq'   : 'CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG',
                    'fss'    : '((((((....((((((((((((....))))))))))))....))))))',
                    # 'cseq'   : 'CTAGGATATGGUUUAUAGGCGGGAACGCCUAUAAACCTAAGTCCTAG', #highest
                    'cseq'   : 'CTAGGATATGGGGGGUUUUUGGGAACAAAAACCCCCCTAAGTCCTAG', #lowest
                    'css'    : '(((((((..((((((((((((....))))))))))))...)))))))'}
        self.options = option.Options(options)

    def _remove_Us(self, seq):
        seq_rna = ""
        for e in seq:
            if e == 'T':
                seq_rna += 'U'
            else:
                seq_rna += e
        return seq_rna

    def _get_motifs_from_seq_and_ss(self, seq, ss):
        parser = secondary_structure_parser.SecondaryStructureParser()
        mg = parser.parse_to_motif_graph(seq, ss)

        start = 0
        motifs = []
        for n in mg.graph.nodes:
            if n.data.mtype == motif_type.TWOWAY and start == 0:
                start = 1
                continue

            if n.data.mtype == motif_type.HAIRPIN:
                break

            if not start:
                continue

            if n.data.mtype == motif_type.HELIX:
                motif = rm.manager.get_bp_step(n.data.end_ids[0])
                motifs.append(motif)
            elif n.data.mtype == motif_type.TWOWAY:
                try:
                    motif = rm.manager.get_motif(end_id=n.data.end_ids[0])
                    motifs.append(motif)
                except exceptions.ResourceManagerException:
                    raise ValueError(
                        "cannot find a motif that corresponds to the sequence: " +
                        motif.sequence() + " and secondary structure: " +
                        motif.dot_bracket() + " for the simulation")
            else:
                raise ValueError(
                    "motif type: " + motif_type.type_to_str(n.data.mtype) + " is not "
                    "supported in tecto simulations.")

        return motifs

    def _get_mset(self):
        fseq = self._remove_Us(self.option('fseq'))
        fss  = self.option('fss')

        cseq = self._remove_Us(self.option('cseq'))
        css  = self.option('css')

        flow_motifs = self._get_motifs_from_seq_and_ss(fseq, fss)
        chip_motifs = self._get_motifs_from_seq_and_ss(cseq, css)

        mt = motif_tree.MotifTree()
        mt.option('sterics', 0)
        m = rm.manager.get_bp_step("GG_LL_CC_RR")
        mt.add_motif(m)
        mt.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")
        mt.add_motif(flow_motifs[1], parent_end_name="A7-A22")

        for i in range(2, len(flow_motifs)):
            mt.add_motif(flow_motifs[i])

        mt.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
        mt.add_motif(chip_motifs[1], parent_end_name="A222-A251")

        for i in range(2, len(chip_motifs)):
            mt.add_motif(chip_motifs[i])

        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt=mt)
        self.mt = mt
        return mset

    def run(self):

        analyse = 0

        ni1 = 2

        self.mset = self._get_mset()
        self.mst = self.mset.to_mst()
        ni2 = self.mst.last_node().index
        # ni2 = 3
        if analyse == 0:
            msec = MotifStateEnsembleConvolution(self.mset,ni1,ni2,grid_size=[21,21,20,21],grid_unit=1.5)
            msec.run()
        else:
            npzf = np.load('run2.npz')
            msec = MotifStateEnsembleConvolution(self.mset,ni1,ni2,grid_size=[9,9,10,9],grid_unit=2.0)
            msec.map.data = npzf['lastdata']
            msec.slidetrack = [npzf['slidetotal']]
            msec.done_ = True
        start_state = self.mst.get_node(ni1).data.cur_state.end_states[0]
        start_matrix = se3.state_to_matrix(start_state)
        print '\nstart_matrix',start_matrix
        end_state = self.mst.get_node(1).data.cur_state.end_states[1]
        end_matrix = se3.state_to_matrix(end_state)
        print '\nend_matrix',end_matrix
        step_matrix = np.dot(la.inv(start_matrix),end_matrix)

        end_matrix_mean = msec.gl.get_mg(ni1,ni2).mean
        print '\n end_matrix_mean',end_matrix_mean
        print '\nstep matrix_ginger',step_matrix
        # end_state = self.mst.get_node(ni2).data.cur_state.end_states[1]
        # end_matrix = se3.state_to_matrix(end_state)
        # step_matrix = np.dot(la.inv(start_matrix), end_matrix)
        # print '\nstep matrix_mst', step_matrix
        pdf = msec.get_prob(step_matrix)
        print 'normalized pdf:', pdf
        if analyse ==1:
            apdf = msec.get_prob_analysis(step_matrix)
            print 'pdf for analysis:',apdf




if __name__ == "__main__":
    args = parse_args()
    opts = vars(args)
    st = SimulateTectos()
    st._get_mset()
    mset = st._get_mset()
    mst = mset.to_mst()
    mt = st.mt
    mt.to_pdb("test.pdb", renumber=1, close_chain=1)
    st.run()



