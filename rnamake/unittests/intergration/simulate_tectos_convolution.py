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

    def run(self):

        # init
        self.map.place_ensemble(self.mset.get_node(self.ni_f_).data)
        self.slidemap(self.ni_f_)

        #iter
        for en in self.mset.tree.nodes:#use iterator, start from one
            print 'currently process node # %d'%en.index
            if en.index<=self.ni_f_:
                continue
            if en.index> self.ni_l_:
                self.done_ = True
                print 'self.done should be updated!!!'
                break
            #if statement
            temp_map = se3.SE3Map(self.grid_size, self.grid_unit)
            temp_map.place_ensemble(en.data)
            self.map = self.map*temp_map
            self.slidemap(en.index)
            del temp_map
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
        return self.map.data[self.map.matrix_to_grid_ndx(m)]

    def slidemap(self,index):
        ni = copy.copy(index)
        meand = self.gl.mgl[ni].mean[:3,3]
        # meand = self.mst_.get_node(ni).data.ref_state.end_states[1].d
        meangrid = np.around(meand/self.grid_unit)
        print 'meand',meand
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
                    'cseq'   : 'CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG',
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
        ni1 = 11
        ni2 = 13
        self.mset = self._get_mset()
        self.mst = self.mset.to_mst()
        msec = MotifStateEnsembleConvolution(self.mset,ni1,ni2,grid_size=[11,11,11,11],grid_unit=3.0)
        msec.run()
        start_state = self.mst.get_node(ni1).data.cur_state.end_states[0]
        start_matrix = se3.state_to_matrix(start_state)
        end_state = self.mst.get_node(ni2).data.cur_state.end_states[1]
        end_matrix = se3.state_to_matrix(end_state)
        step_matrix = np.dot(la.inv(start_matrix),end_matrix)
        pdf = msec.get_prob(step_matrix)
        print pdf



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



