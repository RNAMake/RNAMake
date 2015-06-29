import unittest
import random
import rnamake.motif_state_ensemble_tree as motif_state_ensemble_tree
import rnamake.resource_manager as rm
import rnamake.motif_type as motif_type
import rnamake.util as util
import rnamake.basic_io as basic_io
import rnamake.settings as settings
import rnamake.motif as motif
import rnamake.motif_tree as motif_tree
import rnamake.motif_type as motif_type

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




class MotifStateEnsembleTreeUnittest(unittest.TestCase):

    def test_creation(self):
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree()

    def test_add(self):
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree()
        mse = rm.manager.get_motif_state_ensemble("GG_LL_CC_RR")
        mset.add_ensemble(mse)
        mset.add_ensemble(mse)
        mst = mset.to_mst()
        mt = mst.to_motif_tree()


    def _test_get_mtst(self):
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        me2 = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        me3 = motif_ensemble.MotifEnsemble("GC=AU", 0, 0)
        me4 = motif_ensemble.MotifEnsemble("AU=AU", 0, 0)
        met = motif_ensemble_tree.MotifEnsembleTree()
        met.add_ensemble(me)
        met.add_ensemble(me2)
        met.add_ensemble(me3)
        met.add_ensemble(me4)
        mtst = met.get_mtst()
        print len(mtst.nodes)
        mtst.nodes_to_pdbs()

        h_lib = motif_library.MotifLibrary(motif_type.HELIX)
        mt = motif_tree.MotifTree()
        mt.add_motif(h_lib.get_motif("HELIX.IDEAL.12"), end_index=1)
        mt.write_pdbs("h")

    def _test_get_mtst_2(self):
        pos_bps = ["GC", "CG", "AU", "AU"]
        size = random.randint(10,20)
        bps = []
        for i in range(size):
            bps.append(random.choice(pos_bps))

        names = []
        for i in range(1,size):
            names.append(bps[i-1] + "=" + bps[i])

        met = motif_ensemble_tree.MotifEnsembleTree()
        for n in names:
            me = motif_ensemble.MotifEnsemble(n, 0, 0)
            met.add_ensemble(me)

        mtst = met.get_mtst()
        print len(mtst.nodes)
        mtst.nodes_to_pdbs()

    def _test_mts(self):
        path = settings.RESOURCES_PATH + "prediction/GC=GC.new.me"
        mts_lib = motif_tree_state.MotifTreeStateLibrary(libpath=path, new=1)
        mtst = motif_tree_state.MotifTreeStateTree(sterics=1)

        mts = mts_lib.get_state("GC=GC.2-0")
        mtst.add_state(mts_lib.get_state("GC=GC.2-0"))
        mtst.add_state(mts_lib.get_state("GC=GC.2-0"))
        mtst.add_state(mts_lib.get_state("GC=GC.2-0"))
        mtst.add_state(mts_lib.get_state("GC=GC.2-0"))
        mtst.add_state(mts_lib.get_state("GC=GC.2-0"))
        print len(mtst.nodes)
        mtst.nodes_to_pdbs()
        return
        m = motif.str_to_motif(mts_lib.get_state("GC=GC.4-0").build_string)
        mt = motif_tree.MotifTree()
        mt.add_motif(m, end_index=0, end_flip=0)
        mt.add_motif(m, end_index=0, end_flip=0)
        mt.write_pdbs()
        print len(mt.nodes)

    def _get_chain_pos(self, p, r):
        for i, c in enumerate(p.chains()):
            if r in c.residues:
                return i


    def _problem_one(self):
        twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
        helixs = motif_tree_state.MotifTreeStateLibrary(libpath='HELIX_test.new.me', exclude=['11','10'])
        mtst = motif_tree_state.MotifTreeStateTree()
        #mtst.add_state(helixs.get_state('HELIX.LE.6-0-0-0-0-1-1'))
        mtst.add_state(twoways.get_state('TWOWAY.2VQE.50-0-0-0-0-1-1'))
        node = mtst.add_state(helixs.get_state('HELIX.LE.6-0-0-0-0-1-1'))
        print node
        mtst.nodes_to_pdbs()
        for i, n in enumerate(mtst.nodes):
            basic_io.points_to_pdb('beads.'+str(i)+".pdb", n.beads)


        return mtst

    def _problem_test(self):
        twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
        helixs = motif_tree_state.MotifTreeStateLibrary(libpath='HELIX_test.new.me', exclude=['11','10'])
        mtst = motif_tree_state.MotifTreeStateTree()
        mts_libs = [ helixs, twoways ]

        f = open("test")
        lines = f.readlines()
        f.close()
        pos = 0

        for i, l in enumerate(lines):
            name = l.rstrip()
            if i % 2 == 0 :
                pos = 0
            else:
                pos = 1

            node = mtst.add_state(mts_libs[pos].get_state(name))
            if node is None:
                print "cannot add", name
                exit()
        return mtst




    def _test_mtst_to_met_helix(self):
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


    def _test_mtst_to_met(self):
        #mtst = get_twoway_helix_mts_tree(50)
        #for i, n in enumerate(mtst.nodes):
        #    print i, n.mts.name
        mtst = self._problem_test()
        converter = motif_ensemble_tree.MTSTtoMETConverter()
        mtst2, score = converter.convert(mtst, debug=1)
        print len(mtst2.nodes)

        state_1 = mtst.last_node.active_states()[0]
        state_2 = mtst2.last_node.active_states()[0]
        dist = util.distance(state_1.d, state_2.d)
        print dist

        mtst.nodes_to_pdbs()
        for i, n in enumerate(mtst.nodes):
            basic_io.points_to_pdb('beads.'+str(i)+".pdb", n.beads)


        mtst.to_pdb('test.pdb')
        mtst2.to_pdb('test2.pdb')
        #motif_ensemble_tree.mtst_to_met(mtst)

    def _test_mtst_to_met_exhustive(self):
        return
        twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
        helixs = motif_tree_state.MotifTreeStateLibrary(libpath='HELIX_test.new.me', exclude=['11','10'])

        mtst = motif_tree_state.MotifTreeStateTree()
        for mts1 in helixs.motif_tree_states:
            mtst.add_state(mts1)
            for mts2 in twoways.motif_tree_states:
                node = mtst.add_state(mts2)
                if node is None:
                    continue
                converter = motif_ensemble_tree.MTSTtoMETConverter()
                try:
                    mtst2 = converter.convert(mtst, debug=1)
                except:
                    mtst.remove_node()
                    continue

                state_1 = mtst.last_node.active_states()[0]
                state_2 = mtst2.last_node.active_states()[0]
                dist = util.distance(state_1.d, state_2.d)
                if dist > 5:
                    print dist, mts1.name, mts2.name
                mtst.remove_node()
            mtst.remove_node()






def main():
    unittest.main()

if __name__ == '__main__':
    main()
