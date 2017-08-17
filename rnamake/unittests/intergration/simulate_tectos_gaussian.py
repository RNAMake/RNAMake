from rnamake import base, option, motif_state_ensemble_tree
from rnamake import secondary_structure_parser, motif_type, motif_tree, exceptions
from rnamake import resource_manager as rm
from rnamake import se3util as se3
import numpy as np
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
    def __init__(self, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        #self.mset = self._get_mset()

    def setup_options_and_constraints(self):
        options = { 'fseq'   : 'CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG',
                    'fss'    : '((((((....((((((((((((....))))))))))))....))))))',
                    # 'cseq'   : 'CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG',
                    'cseq'   : 'CTAGGATATGGUUUAUAGGCGGGAACGCCUAUAAACCTAAGTCCTAG',
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
        self.mset = self._get_mset()
        self.mst = self.mset.to_mst()

        ni1 = 2
        # ni2 = self.mst.last_node().index
        ni2 = 3
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

        '''testing chi'''
        from scipy import linalg as la
        target_chi = se3.matrix_to_chi(np.dot(la.inv(mg.mean),target_matrix))
        print target_chi
        print mg.eval(target_chi)
        test_chi = target_chi
        # test_chi = np.array([0,0,0,0,0,0])
        print 'PDF at destination: ', mg.eval(test_chi)
        test_chi[3] += 2
        print 'PDF deviated x=2 A: ', mg.eval(test_chi)
        test_chi[3] -=2
        test_chi[4] += 2
        print 'PDF deviated y=2: ', mg.eval(test_chi)
        test_chi[4] -= 2
        test_chi[5] += 2
        print 'PDF deviated z=2: ', mg.eval(test_chi)
        test_chi[5] -= 2
        # test_chi[4] +=2
        # test_chi = np.array([0,0,0,0,0,0]).astype('float')
        test_n = 100
        test_grid = np.mgrid[-0.2:0.2:test_n*1j]
        from matplotlib import pylab as plt
        plt.figure('PDF over v_x')
        v = np.zeros([test_n])
        for x in range(test_n):
            temp_chi = test_chi.copy()
            temp_chi[0]+=test_grid[x]
            v[x]=mg.eval(temp_chi)
        plt.plot(test_grid,v)

        plt.figure('PDF over v_y')
        v = np.zeros([test_n])
        for x in range(test_n):
            temp_chi = test_chi.copy()
            temp_chi[1] += test_grid[x]
            v[x] = mg.eval(temp_chi)
        plt.plot(test_grid,v)

        plt.figure('PDF over v_z')
        v = np.zeros([test_n])
        for x in range(test_n):
            temp_chi = test_chi.copy()
            temp_chi[2] += test_grid[x]
            v[x] = mg.eval(temp_chi)
        plt.plot(test_grid,v)

        test_grid = np.mgrid[-2:2:test_n * 1j]
        plt.figure('PDF over x')
        v = np.zeros([test_n])
        for x in range(test_n):
            temp_chi = test_chi.copy()
            temp_chi[3] += test_grid[x]
            v[x] = mg.eval(temp_chi)
        plt.plot(test_grid, v)

        plt.figure('PDF over y')
        v = np.zeros([test_n])
        for x in range(test_n):
            temp_chi = test_chi.copy()
            temp_chi[4] += test_grid[x]
            v[x] = mg.eval(temp_chi)
        plt.plot(test_grid, v)

        plt.figure('PDF over z')
        v = np.zeros([test_n])
        for x in range(test_n):
            temp_chi = test_chi.copy()
            temp_chi[5] += test_grid[x]
            v[x] = mg.eval(temp_chi)
        plt.plot(test_grid, v)

        plt.show()

        '''integrated probability'''
        vx,vy,vz,x,y,z = target_chi[0],target_chi[1],target_chi[2],\
        target_chi[3],target_chi[4],target_chi[5]
        from scipy import integrate as jf
        def wrapper(x0, x1, x2):
            return mg.eval(np.array([x0, x1, x2, x,y,z]))
        prob = jf.nquad(wrapper,
                 [[vx-0.2,vx+0.2],[vy-0.2,vy+0.2],[vz-0.2,vz+0.2]],opts={'epsrel':0.2})
        print 'result',prob
        factor = prob[0]/mg.eval(target_chi)
        def wrapper(x0,x1,x2):
            return mg.eval(np.array([vx,vy,vz,x0,x1,x2]))
        prob = jf.nquad(wrapper,
                 [[x - 2, x + 2], [y - 2, y + 2], [z - 2, z + 2]], opts={'epsabs': 0.1, 'epsrel': 0.2})
        print prob
        print 'further result',prob[0]*factor
        # print mg.eval(target_chi)*32.*8*np.pi**2


        '''comparison of mean'''
        print 'ni2 resultant mean\n',mg.mean,'\n'
        state_node = self.mst.get_node(ni2).data
        print 'ni2 resultant state',state_node.get_end_state(state_node.end_name(1))
        # for x in range(0,ni2+1):
        print 'the state to align to ', state_node.get_end_state(state_node.end_name(1))
        for x in range(0,ni2+1):
            state_node = self.mst.get_node(x).data
            print 'node %d end 0 current\n'%x, state_node.get_end_state(state_node.end_name(1))
        for x in range(0,ni2+1):
            state_node = self.mst.get_node(x).data
            print 'node %d end 0 ref\n'%x,state_node.ref_state.end_states[1]







if __name__ == "__main__":
    args = parse_args()
    opts = vars(args)
    st = SimulateTectos()
    st._get_mset()
    mt = st.mt
    mt.to_pdb("test.pdb", renumber=1, close_chain=1)
    st.run()