import unittest
import random
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree
import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.motif_tree_state as motif_tree_state
import rnamake.motif_type as motif_type

def get_twoway_helix_mts_tree(size=2):
    twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
    helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
    me_libs = [helixs, twoways]
    mtst = motif_tree_state.MotifTreeStateTree()
    pos = 0
    i = 0
    count = 0
    while i < size:
        if i % 2 == 0:
            pos = 0
        else:
            pos = 1
        mts = random.choice(me_libs[pos].motif_tree_states)
        node = mtst.add_state(mts)
        if node is not None:
            i += 1
        count += 1
        if count > 1000:
            break
    return mtst




class MotifEnsembleTreeUnittest(unittest.TestCase):

    def test_creation(self):
        met = motif_ensemble_tree.MotifEnsembleTree()

    def test_add(self):
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        met = motif_ensemble_tree.MotifEnsembleTree()
        met.add_ensemble(me)

    def test_get_mtst(self):
        return
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        me2 = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        met = motif_ensemble_tree.MotifEnsembleTree()
        met.add_ensemble(me)
        met.add_ensemble(me2)
        mtst = met.get_mtst()
        print len(mtst.nodes)
        mtst.nodes_to_pdbs()

    def get_chain_pos(self, p, r):
        for i, c in enumerate(p.chains()):
            if r in c.residues:
                return i


    def _problem_one(self):
        twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
        helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
        mtst = motif_tree_state.MotifTreeStateTree()
        mtst.add_state(helixs.get_state('HELIX.LE.14-0-0-1-14-0-0'))
        mtst.add_state(twoways.get_state('TWOWAY.2VQE.13-0-0-0-0-1-1'))
        return mtst


    def test_mtst_to_met(self):
        #mtst = get_twoway_helix_mts_tree(2)
        mtst = self._problem_one()
        motif_ensemble_tree.mtst_to_met(mtst)

        exit()

        p = mtst.to_pose()
        p.to_pdb('test.pdb')
        n1 = p.nodes[1]
        chain_ends = []
        for c in n1.motif.chains():
           chain_ends.append(c.first())
           if len(c) > 1:
               chain_ends.append(c.last())

        start_chain = None
        start_res = None
        chain_pos = 0
        for r in chain_ends:
            r_new = p.get_residue(uuid=r.uuid)
            if r_new is None:
                continue

            found = 0
            for i, c in enumerate(p.chains()):
                if r_new == c.first():
                    found = 1
                    start_res = r_new
                    start_chain = c
                    chain_pos = i
                    break

        designed_sequence = p.optimized_sequence()
        #print p.secondary_structure(
        print
        print designed_sequence
        print p.secondary_structure()

        bps_str = []
        seen = []
        last_res_i = 0
        for i, r in enumerate(start_chain.residues):
            last_res_i = i
            r_new = p.nodes[1].motif.get_residue(uuid=r.uuid)
            if r_new is None:
                break
            #print r.num, chain_pos
            res1 = designed_sequence[r.num-1 + chain_pos]
            bps = p.get_basepair(uuid1=r.uuid)
            for bp in bps:
                if bp.res1 == r:
                    res2 =  designed_sequence[bp.res2.num-1 + \
                                              self.get_chain_pos(p, bp.res2)]
                else:
                    res2 =  designed_sequence[bp.res1.num-1 + \
                                              self.get_chain_pos(p, bp.res1)]

            bps_str.append(res1+res2)

        steps = []
        for i in range(1, len(bps_str)):
            steps.append(bps_str[i-1]+"="+bps_str[i])

        print steps
        met = motif_ensemble_tree.MotifEnsembleTree()
        for s in steps:
            me = motif_ensemble.MotifEnsemble(s, 0, 0)
            met.add_ensemble(me)

        mtst2 = met.get_mtst()
        mtst2.to_pdb("test2.pdb")






def main():
    unittest.main()

if __name__ == '__main__':
    main()
