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
                    # 'cseq'   : cseq,
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
        # self.mset = self._get_mset()
        self.mst = self.mset.to_mst()

        ni1 = 2
        ni2 = self.mst.last_node().index
        # ni2 = 10
        # default ei = 1
        mgl = se3.MotifGaussianList(self.mset)
        mg = mgl.get_mg(ni1,ni2)
        state_node = self.mst.get_node(ni1).data
        mg.mean = np.dot(se3.state_to_matrix(
            state_node.get_end_state(state_node.end_name(0))
        ),mg.mean)
        state_node =self.mst.get_node(1).data
        target_matrix = se3.state_to_matrix(state_node.get_end_state(state_node.end_name(1)))

        # '''printing out the topology'''
        # print self.mst.to_pretty_str()
        '''matrix for double gaussian'''
        # cutoff = 4.7
        # radius = 10.0
        # radian = cutoff/radius
        # ele_v = 2.0/radian**2
        # ele_x = 2.0/cutoff**2
        # TAU = np.diag([ele_v*2,ele_v,ele_v/5,ele_x,ele_x,ele_x])
        cutoff = 5
        radius = 10.0
        radian = cutoff/radius
        ele_v = 2.0/radian**2
        ele_x = 2.0/cutoff**2
        TAU = np.diag([ele_v*2,ele_v,ele_v,ele_x,ele_x,ele_x])
        print 'TAU diag'
        print [ele_v*2,ele_v,ele_v,ele_x,ele_x,ele_x]
        '''testing chi'''
        from scipy import linalg as la
        target_chi = se3.matrix_to_chi(np.dot(la.inv(mg.mean),target_matrix))
        print target_chi
        print 'PDF at target: ',mg.eval(target_chi)
        res = mg.eval_double(target_chi,TAU)
        print 'P by double Gaussian: ',res
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
    def __init__(self):
        #
        # self.datajoe = np.genfromtxt('/home/zhuoyu/Downloads/exhustive_helices.results',
        #                         skip_header=1,usecols=(0,1),dtype=['S128','float'])
        # self.data_count = self.datajoe.shape[0]

        self.dataexp = np.genfromtxt('/home/zhuoyu/Downloads/helical_variation.csv',skip_header=1,
                                     usecols=(0,1),dtype=['S128','float'],delimiter=',')
        self.data_count = self.dataexp.shape[0]

        self.datajoe = np.genfromtxt('/home/zhuoyu/Downloads/exhustive_helices_TU.results',skip_header=1,
                                     usecols=(0,1),dtype=['S128','float'])
        self.wtctjoe = 13060.0

        self.n=self.data_count
        self.cmpdata = np.zeros([3,self.n])
        self.cmpdata[:] = np.nan
        self.k1 =0
        wtst = SimulateTectos(cseq='CTAGGATATGGGGUAGGUGCGGGAACGCACCUACCCCTAAGTCCTAG')
        self.wtprob = wtst.run()
        self.wtdg = -10.583
    def ud(self,x):
        # i = np.random.randint(self.data_count)
        i = x
        cseq,dg = self.dataexp[i]
        st = SimulateTectos(cseq=str(cseq))

        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)
        # st._get_mset()
        # mt = st.mt
        # mt.to_pdb("test.pdb", renumber=1, close_chain=1)
        ndxjoe = np.where(self.datajoe['f0'] == cseq)[0]
        if len(ndxjoe) == 1:
            probjoe = self.datajoe[ndxjoe]['f1'][0]
        # print self.datajoe[np.where(self.datajoe['f0'] == cseq)]['f0'][0]
        # print cseq
        # print probjoe
            ddgjoe  = -kBT*6.02/4.1868/100*np.log(probjoe/self.wtctjoe)
            dgjoe = ddgjoe + self.wtdg
        else:
            dgjoe = np.nan
            print 'data unmatched in Joe\'s and exp'
        simulprob = st.run()
        simulddg = -kBT*6.02/4.1868/100*np.log(simulprob/self.wtprob)
        simuldg = simulddg+self.wtdg
        self.cmpdata[0,self.k1],self.cmpdata[1,self.k1],self.cmpdata[2,self.k1] =simuldg,dg,dgjoe
        self.k1 +=1


class RunMeVarLen():
    def __init__(self):
        self.dataexp = np.genfromtxt('/home/zhuoyu/Downloads/all_lengths.csv', skip_header=1,
                                usecols=(0, 1,2,3,4), dtype=['S128','S128','S128','S128','float'],
                                     delimiter=',')
        self.data_count = self.dataexp.shape[0]

        self.datajoe = np.genfromtxt('/home/zhuoyu/Downloads/org_wcall_results.csv', skip_header=1,
                                     usecols=(1,3,8), dtype=['S128', 'S128','float'],delimiter=',')
        self.wtctjoe = 13060.0


        self.n = self.data_count
        self.cmpdata=np.zeros([3,self.n])
        self.cmpdata[:]=np.nan
        self.k1 = 0

        wtst = SimulateTectos(cseq='CTAGGATATGGGGUAGGUGCGGGAACGCACCUACCCCTAAGTCCTAG')
        self.wtprob = wtst.run()
        self.wtdg = -10.583

    def ud(self,x,query_flen,query_clen):
        i = x
        fseq,fss,cseq,css,dg = self.dataexp[i]
        flen = 10+(len(str(fss))-len("CUAGGAAUCUGGAAGUACCGAGGAAACUCGGUACUUCCUGUGUCCUAG"))/2
        clen = 10+(len(str(css))-len("CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG"))/2
        if flen != query_flen or clen != query_clen:
            return
        # else:
        #     self.cmpdata[0, self.k1], self.cmpdata[1, self.k1], self.cmpdata[2, self.k1] = i-1,i,i+1
        #     self.k1 += 1
        #     return
        st = SimulateTectos(fseq=str(fseq),fss=str(fss),cseq=str(cseq),css=str(css))
        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)
        # st._get_mset()
        # mt = st.mt
        # mt.to_pdb("test.pdb", renumber=1, close_chain=1)
        ndxjoe = np.where(np.logical_and(self.datajoe['f0'] == fseq,self.datajoe['f1']==cseq))[0]
        if len(ndxjoe) == 1:
            ddgjoe = self.datajoe[ndxjoe]['f2'][0]
            # print self.datajoe[np.where(self.datajoe['f0'] == cseq)]['f0'][0]
            # print cseq
            # print probjoe
            # ddgjoe = -kBT * 6.02 / 4.1868 / 100 * np.log(probjoe / self.wtctjoe)
            dgjoe = ddgjoe + self.wtdg
        else:
            dgjoe = np.nan
            print 'data unmatched in Joe\'s and exp'
        simulprob = st.run()
        simulddg = -kBT * 6.02 / 4.1868 / 100 * np.log(simulprob / self.wtprob)
        simuldg = simulddg + self.wtdg
        self.cmpdata[0, self.k1], self.cmpdata[1, self.k1], self.cmpdata[2, self.k1] = simuldg, dg, dgjoe
        self.k1 += 1





if __name__ == "__main__":
    # args = parse_args()
    # opts = vars(args)
    runme =  RunMe()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # lenlist = [(10,11)]
    lenlist = [(10,9),(10,10),(10,11),(11,9),(11,10),(9,9),(9,10)]
    # colorlist = cm.rainbow(np.mgrid[0:1:2*len(lenlist)*1j])
    sctlist_mine = []
    sctlist_joe =[]
    len_split = []
    # for ndx_len in enumerate(lenlist):

    # for i,ndx_len in enumerate(lenlist):
    #     len_split.append(runme.k1)
    #     k2 = runme.k1
        # runme.cmpdata[:]=np.nan
        # for x in range(runme.n):
        #     runme.ud(x,*ndx_len)
        # sct_mine= plt.scatter(runme.cmpdata[1, k2:], runme.cmpdata[0, k2:],color=colorlist[i])
        # sct_joe=plt.scatter(runme.cmpdata[1, k2:], runme.cmpdata[2, k2:], marker='^',color=colorlist[i+len(lenlist)])
        # sctlist_mine.append(sct_mine)
        # sctlist_joe.append(sct_joe)
    # for x in range(runme.n):
    #     runme.ud(x)
    npzf = np.load('gaussian_1010.npz')
    runme.cmpdata = npzf['data']
    sct_mine= plt.scatter(runme.cmpdata[1, :], runme.cmpdata[0, :])
    plt.plot(np.mgrid[-11:-10:100j], np.mgrid[-11:-10:100j], 'k')
    plt.figure()
    sct_joe=plt.scatter(runme.cmpdata[1, :], runme.cmpdata[2, :], marker='^')
    # np.savez('gaussian_1010',data=runme.cmpdata)
    # legend1 =plt.legend(tuple(sctlist_joe),('flow %d chip %d' % x for x in lenlist),loc='lower right',
    #                     scatterpoints=1,title='MC')
    # legend2=plt.legend(tuple(sctlist_mine), ('flow %d chip %d' % x for x in lenlist), scatterpoints=1,
    #            loc='upper left',title='Gaussian')
    # legend1.draggable()
    # legend2.draggable()
    # plt.gca().add_artist(legend1)
    plt.plot(np.mgrid[-11:-10:100j],np.mgrid[-11:-10ls
    len_split:100j],'k')
    ax.set_xlabel('Exp $\Delta G$(kcal/mol)')
    ax.set_ylabel('Predicted $\Delta G$(kcal/mol)')
    plt.show()





