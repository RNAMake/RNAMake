from rnamake import base, option, motif_state_ensemble_tree
from rnamake import secondary_structure_parser, motif_type, motif_tree, exceptions
from rnamake import resource_manager as rm
from rnamake import se3util as se3
import numpy as np
from matplotlib import pylab as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
# import warnings
# warnings.filterwarnings('error')


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
    def __init__(self,**options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.mset = self._get_mset()

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
        self.mst = self.mset.to_mst()
        # this is only valid for tectoRNAs, of which the first index is 2, and the last index is the motif to be
        # aligned to the first tetraloop. The first tetraloop has the index 1, and the last motif should be aligned
        # to the end[1] of it.
        ni1 = 2
        # first index
        ni2 = self.mst.last_node().index
        # last index

        # default ei = 1
        mgl = se3.MotifGaussianList(self.mset)
        # construct a MotifGaussianList from the mset
        mg = mgl.get_mg(ni1,ni2)
        # get the resultant MotifGaussian from ni1 to ni2
        state_node = self.mst.get_node(ni1).data
        mg.mean = np.dot(se3.state_to_matrix(
            state_node.get_end_state(state_node.end_name(0))
        ),mg.mean)
        # align the former MotifGaussian to the state of the starting basepair in the MotifTree, so that all the
        # coordinates are w.r.t. the origin of the MotifTree, expected to be the end[0] of the motif[0]
        state_node =self.mst.get_node(1).data
        target_matrix = se3.state_to_matrix(state_node.get_end_state(state_node.end_name(1)))
        # get the state matrix of the basepair for the last motif to be aligned to, namely the 'target'

        # '''printing out the topology'''
        # print self.mst.to_pretty_str()
        '''matrix for double gaussian'''
        # below is the parameters used to construct the window Gaussian.

        # This is a set of parameters that work decently and are simple to understand
        # cutoff = 4.7
        # radius = 10.0
        # radian = cutoff/radius
        # ele_v = 2.0/radian**2
        # ele_x = 2.0/cutoff**2
        # TAU = np.diag([ele_v*2,ele_v,ele_v/5,ele_x,ele_x,ele_x])

        # This is fine-tuned set of parameters
        cutoff = 4.5
        # The cutoff radius in Angstrom
        radius = 10.0
        # The radius of the basepair, in Angstrom
        radian = cutoff/radius
        ele_v = 2.0/radian**2
        # The elements of the TAU matrix is calculated based on that the 'cutoff' of a Gaussian is defined
        # by its decaying to e^(-1) to the peak.
        cutoffv = np.array([4.5, 4.5, 4.0])
        # the cutoff radius specified in each direction
        ele_xv = 2.0 / cutoffv ** 2

        TAU = np.diag([ele_v, ele_v, ele_v, ele_xv[0], ele_xv[1], ele_xv[2]])

        '''testing chi'''
        from scipy import linalg as la
        target_chi = se3.matrix_to_chi(np.dot(la.inv(mg.mean), target_matrix)) + \
                     np.array([0.0, 0.0, 0.0, 0 * np.sqrt(2.0 / TAU[3, 3]), 0.0, 2.0 * np.sqrt(2.0 / TAU[5, 5])])
        # we may slightly adjust that target_chi, which is natrual, since we may add an exponential factor before the
        # Gaussian, representing the repulsion of the other parts of the RNA, and the result would always be a Gaussian
        # with shifted mean.
        print target_chi
        print 'node1 end2 chi', se3.matrix_to_chi(
            np.dot(la.inv(mg.mean), se3.state_to_matrix(self.mst.get_node(1).data.cur_state.end_states[2])))
        print 'PDF at target: ',mg.eval(target_chi)
        res = mg.eval_double(target_chi,TAU)
        print 'P by double Gaussian: ',res

        '''Below is only used for all kinds of debugging, and plotting'''
        #
        # test_chi = target_chi
        # # test_chi = np.array([0,0,0,0,0,0])
        # print 'PDF at destination: ', mg.eval(test_chi)
        # test_chi[3] += 2
        # print 'PDF deviated x=2 A: ', mg.eval(test_chi)
        # test_chi[3] -=2
        # test_chi[4] += 2
        # print 'PDF deviated y=2: ', mg.eval(test_chi)
        # test_chi[4] -= 2
        # test_chi[5] += 2
        # print 'PDF deviated z=2: ', mg.eval(test_chi)
        # test_chi[5] -= 2
        # # test_chi[4] +=2
        # # test_chi = np.array([0,0,0,0,0,0]).astype('float')
        # test_n = 100
        # test_grid = np.mgrid[-0.2:0.2:test_n*1j]
        #
        # plt.figure('PDF over v_x')
        # v = np.zeros([test_n])
        # for x in range(test_n):
        #     temp_chi = test_chi.copy()
        #     temp_chi[0]+=test_grid[x]
        #     v[x]=mg.eval(temp_chi)
        # plt.plot(test_grid,v)
        #
        # plt.figure('PDF over v_y')
        # v = np.zeros([test_n])
        # for x in range(test_n):
        #     temp_chi = test_chi.copy()
        #     temp_chi[1] += test_grid[x]
        #     v[x] = mg.eval(temp_chi)
        # plt.plot(test_grid,v)
        #
        # plt.figure('PDF over v_z')
        # v = np.zeros([test_n])
        # for x in range(test_n):
        #     temp_chi = test_chi.copy()
        #     temp_chi[2] += test_grid[x]
        #     v[x] = mg.eval(temp_chi)
        # plt.plot(test_grid,v)
        #
        # test_grid = np.mgrid[-5:5:test_n * 1j]
        # plt.figure('PDF over x')
        # v = np.zeros([test_n])
        # for x in range(test_n):
        #     temp_chi = test_chi.copy()
        #     temp_chi[3] += test_grid[x]
        #     v[x] = mg.eval(temp_chi)
        # plt.plot(test_grid, v)
        #
        # plt.figure('PDF over y')
        # v = np.zeros([test_n])
        # for x in range(test_n):
        #     temp_chi = test_chi.copy()
        #     temp_chi[4] += test_grid[x]
        #     v[x] = mg.eval(temp_chi)
        # plt.plot(test_grid, v)
        #
        # plt.figure('PDF over z')
        # v = np.zeros([test_n])
        # for x in range(test_n):
        #     temp_chi = test_chi.copy()
        #     temp_chi[5] += test_grid[x]
        #     v[x] = mg.eval(temp_chi)
        # plt.plot(test_grid, v)
        #
        # plt.show()
        #
        # '''integrated probability'''
        # vx,vy,vz,x,y,z = target_chi[0],target_chi[1],target_chi[2],\
        # target_chi[3],target_chi[4],target_chi[5]
        # from scipy import integrate as jf
        # def wrapper(x0, x1, x2):
        #     return mg.eval(np.array([x0, x1, x2, x,y,z]))
        # prob = jf.nquad(wrapper,
        #          [[vx-0.2,vx+0.2],[vy-0.2,vy+0.2],[vz-0.2,vz+0.2]],opts={'epsrel':0.2})
        # print 'result',prob
        # factor = prob[0]/mg.eval(target_chi)
        # def wrapper(x0,x1,x2):
        #     return mg.eval(np.array([vx,vy,vz,x0,x1,x2]))
        # prob = jf.nquad(wrapper,
        #          [[x - 2, x + 2], [y - 2, y + 2], [z - 2, z + 2]], opts={'epsabs': 0.1, 'epsrel': 0.2})
        # print prob
        # print 'further result',prob[0]*factor
        # # print mg.eval(target_chi)*32.*8*np.pi**2
        #
        #
        # '''comparison of mean'''
        # print 'ni2 resultant mean\n',mg.mean,'\n'
        # state_node = self.mst.get_node(ni2).data
        # print 'ni2 resultant state',state_node.get_end_state(state_node.end_name(1))
        # # for x in range(0,ni2+1):
        # # print 'the state to align to ', state_node.get_end_state(state_node.end_name(1))
        # for x in range(0,ni2):
        #     state_node = self.mst.get_node(x).data
        #     print 'node %d name'%x,state_node.name()
        #     print 'node %d end 0 current\n'%x, state_node.get_end_state(state_node.end_name(0))
        # state_node = self.mst.get_node(2).data
        # node2state = state_node.get_end_state(state_node.end_name(0))
        # print 'target node w.r.t node 2'
        # for x in range(2,ni2+1):
        #     print 'node %d mean w.r.t node2 \n'%x,mgl.get_mg(2,x).mean

        # for x in range(0,ni2+1):
        #     state_node = self.mst.get_node(x).data
        #     print 'node %d end 0 ref\n'%x,state_node.ref_state.end_states[0]
        return res




class RunMe():
    """
    The object to run the calculation with only one length
    """
    def __init__(self):

        self.dataexp = np.genfromtxt('/home/zhuoyu/Downloads/helical_variation.csv',skip_header=1,
                                     usecols=(0,1),dtype=['S128','float'],delimiter=',')
        self.data_count = self.dataexp.shape[0]
        # open experimental data file

        self.datajoe = np.genfromtxt('/home/zhuoyu/Downloads/exhustive_helices_TU.results',skip_header=1,
                                     usecols=(0,1),dtype=['S128','float'])
        self.wtctjoe = 13060.0
        # open MC data file

        self.n=self.data_count

        self.cmpdata = np.zeros([3,self.n])
        self.cmpdata[:] = np.nan
        # initialize an array to hold the data to be compared

        self.k1 =0
        # initialize the index of loop

        wtst = SimulateTectos(cseq='CTAGGATATGGGGUAGGUGCGGGAACGCACCUACCCCTAAGTCCTAG')
        # get the wild type simulation
        self.wtprob = wtst.run()
        # get the wild type probability
        self.wtdg = -10.583
        # get the wild type free energy
    def ud(self,x):
        # i = np.random.randint(self.data_count)
        i = x
        cseq,dg = self.dataexp[i]
        # get the experimental chip seq and dg
        st = SimulateTectos(cseq=str(cseq))
        # set up the object with the chip seq to run

        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

        ndxjoe = np.where(self.datajoe['f0'] == cseq)[0]
        # get the index of the chipseq in the array of MC data
        if len(ndxjoe) == 1:
            # if there is exactly one matching sequence, as it should be
            probjoe = self.datajoe[ndxjoe]['f1'][0]
            # get the probability from the MC data
            ddgjoe  = -kBT*6.02/4.1868/100*np.log(probjoe/self.wtctjoe)
            # get the ddg from the MC data
            dgjoe = ddgjoe + self.wtdg
            # calculate the dg from ddg
        else:
            # if there is multiple, or no matching sequence, there must be some problem. there then is set a
            # placeholder of NaN to suppress plotting it
            dgjoe = np.nan
            print 'data unmatched in Joe\'s and exp'

        simulprob = st.run()
        # the probability from Gaussian is got
        simulddg = -kBT*6.02/4.1868/100*np.log(simulprob/self.wtprob)
        # the ddg is calculated from probabilty
        simuldg = simulddg+self.wtdg
        # then the dg
        self.cmpdata[0,self.k1],self.cmpdata[1,self.k1],self.cmpdata[2,self.k1] =simuldg,dg,dgjoe
        # append all the dg's in the array, for plotting and analysis
        self.k1 +=1
        # increment the loop index


class RunMeVarLen():
    """
    This is the class for Gaussian calculation of variable lengths of sequences
    """
    def __init__(self):
        self.dataexp = np.genfromtxt('/home/zhuoyu/Downloads/all_lengths.csv', skip_header=1,
                                usecols=(0, 1,2,3,4), dtype=['S128','S128','S128','S128','float'],
                                     delimiter=',')
        self.data_count = self.dataexp.shape[0]
        # load the experimental data

        self.datajoe = np.genfromtxt('/home/zhuoyu/Downloads/org_wcall_results.csv', skip_header=1,
                                     usecols=(1,3,8), dtype=['S128', 'S128','float'],delimiter=',')
        # load the MC data
        self.wtctjoe = 13060.0
        # specify the MC count of wildtype


        self.n = self.data_count
        self.cmpdata=np.zeros([3,self.n])
        self.cmpdata[:]=np.nan
        # initialize the array for data to be compared
        self.k1 = 0
        # initialize the loop index

        wtst = SimulateTectos(cseq='CTAGGATATGGGGUAGGUGCGGGAACGCACCUACCCCTAAGTCCTAG')
        self.wtprob = wtst.run()
        # get the wild type probability
        self.wtdg = -10.583
        # specify the wild type dg

    def ud(self,x,query_flen,query_clen):
        i = x
        fseq,fss,cseq,css,dg = self.dataexp[i]
        # get the sequences and secondary structures and dg's to be calculated
        flen = 10+(len(str(fss))-len("CUAGGAAUCUGGAAGUACCGAGGAAACUCGGUACUUCCUGUGUCCUAG"))/2
        clen = 10+(len(str(css))-len("CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG"))/2
        # get the flow piece and chip piece lengths
        if flen != query_flen or clen != query_clen:
            # if the flow piece and chip piece lengths are not the length to be queried this time, just skip.
            # TODO this kind of filtering flow piece and chip piece lengths is inefficient, and may cause problem
            # when the dataset is huge (although unlikely, since the dataset is generated from experiment, and
            # experiment should be much slower than computer)
            self.fin = 0
            # set the finished flag to 0(no calculation run this time)
            return
        st = SimulateTectos(fseq=str(fseq),fss=str(fss),cseq=str(cseq),css=str(css))
        # construct the simulation object
        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)
        ndxjoe = np.where(np.logical_and(self.datajoe['f0'] == fseq,self.datajoe['f1']==cseq))[0]
        # query the index of the sequence in MC data
        if len(ndxjoe) == 1:
            # There should be exactly 1 data matching
            ddgjoe = self.datajoe[ndxjoe]['f2'][0]
            dgjoe = ddgjoe + self.wtdg
        else:
            dgjoe = np.nan
            print 'data unmatched in Joe\'s and exp'
        simulprob = st.run()
        simulddg = -kBT * 6.02 / 4.1868 / 100 * np.log(simulprob / self.wtprob)
        simuldg = simulddg + self.wtdg
        # get the Gaussain dg
        self.cmpdata[0, self.k1], self.cmpdata[1, self.k1], self.cmpdata[2, self.k1] = simuldg, dg, dgjoe
        # append the Gaussian, experimental, MC dg to the array
        self.k1 += 1
        # increment loop index
        self.fin = 1
        # set finish flag to 1






if __name__ == "__main__":
    runme = RunMeVarLen()
    # if to run under the single length dataset, alter this to RunMe() and modify accordingly the content below
    lenlist = [(10,9),(10,10),(10,11),(11,9),(11,10),(9,9),(9,10)]
    # list of lengths to run
    colorlist = cm.rainbow(np.mgrid[0:1:len(lenlist) * 1j])
    # colors for different lengths
    sctlist_mine = []
    # list of scatter plots of Gaussian
    sctlist_joe =[]
    # list of scatter plots of MC
    len_split = []
    # list of indices splitting data for different length, or the indices starting data for each length.
    # For example, in cmpdata, indices 0 thru 9 stores length (10,9), the first element in the list would be 0 , and
    # the second 10.
    # initialize the plots
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    for i, ndx_len in enumerate(lenlist):
        len_split.append(runme.k1)
        k2 = runme.k1
        # k2 is for selectively plotting data of a certain length
        runme.cmpdata[:] = np.nan
        for x in range(runme.n):
            runme.ud(x, *ndx_len)
            if runme.fin == 1:
                break
        sct_mine = ax1.scatter(runme.cmpdata[1, k2:], runme.cmpdata[0, k2:], color=colorlist[i])
        # make the Gaussian plot
        sct_joe = ax2.scatter(runme.cmpdata[1, k2:], runme.cmpdata[2, k2:], marker='^', color=colorlist[i])
        # make the plot of MC data for comparison
        sctlist_mine.append(sct_mine)
        sctlist_joe.append(sct_joe)

    ax1.plot(np.mgrid[-11:-7:100j], np.mgrid[-11:-7:100j], 'k')
    # plot the reference line of y=x

    legend1 = ax2.legend(tuple(sctlist_joe), ('flow %d chip %d' % x for x in lenlist), loc='lower right',
                         scatterpoints=1, title='MC')
    legend2 = ax1.legend(tuple(sctlist_mine), ('flow %d chip %d' % x for x in lenlist), scatterpoints=1,
                         loc='upper left', title='Gaussian')
    # add the legends

    legend1.draggable()
    legend2.draggable()
    # let the legends float

    ax2.plot(np.mgrid[-11:-7:100j], np.mgrid[-11:-7:100j], 'k')
    # Ah, another reference line

    ax1.set_xlabel('Exp $\Delta G$(kcal/mol)')
    ax2.set_xlabel('Exp $\Delta G$(kcal/mol)')
    ax1.set_ylabel('Predicted $\Delta G$(kcal/mol)')
    ax2.set_ylabel('Predicted $\Delta G$(kcal/mol)')
    # Set the labels

    plt.show()
    # Launch the plot.
