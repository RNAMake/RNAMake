import unittest
import random
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree
import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.motif_tree_state as motif_tree_state
import rnamake.motif_type as motif_type
import rnamake.util as util

def get_twoway_helix_mts_tree(size=2):
    twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
    helixs = motif_tree_state.MotifTreeStateLibrary(libpath='HELIX_test.new.me', exclude=['11','10'])
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
        mtst.add_state(helixs.get_state('HELIX.LE.10-0-0-0-0-1-1'))
        mtst.add_state(twoways.get_state('TWOWAY.2VQE.26-0-0-0-0-1-0'))
        #mtst.add_state(helixs.get_state('HELIX.LE.18-0-0-0-0-1-0'))
        return mtst


    def test_mtst_to_met_helix(self):
        return
        helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX, exclude=['11','10'])

        for mts in helixs.motif_tree_states:
            mtst = motif_tree_state.MotifTreeStateTree()
            mtst.add_state(mts)

            converter = motif_ensemble_tree.MTSTtoMETConverter()
            mtst2 = converter.convert(mtst, debug=1)

            d1 =mtst.last_node.active_states()[0].d
            d2 = mtst2.last_node.active_states()[0].d

            if mts.name == 'HELIX.LE.8-0-0-0-0-1-0':
                mtst.to_pdb('test.pdb')
                mtst2.to_pdb('test2.pdb')
                exit()

            dist = util.distance(d1, d2)
            if dist > 2:
                print mts.name, dist


    def test_mtst_to_met(self):
        #TWOWAY.1JBS.0-0-0-0-0-1-1
        mtst = get_twoway_helix_mts_tree(2)
        for n in mtst.nodes:
            print n.mts.name
        #mtst = self._problem_one()
        converter = motif_ensemble_tree.MTSTtoMETConverter()
        mtst2 = converter.convert(mtst, debug=1)

        mtst.to_pdb('test.pdb')
        mtst2.to_pdb('test2.pdb')
        #motif_ensemble_tree.mtst_to_met(mtst)





def main():
    unittest.main()

if __name__ == '__main__':
    main()
