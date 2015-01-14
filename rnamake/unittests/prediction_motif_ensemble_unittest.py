import unittest
import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree
import rnamake.motif_tree_state as motif_tree_state
import rnamake.motif_tree as motif_tree
import rnamake.settings as settings
import rnamake.util as util
import rnamake.basic_io as basic_io
import rnamake.motif as motif
import rnamake.secondary_structure_tree as ss_tree

class MotifEnsembleUnittest(unittest.TestCase):

    def test_creation(self):
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        me2 = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        met = motif_ensemble_tree.MotifEnsembleTree(me)
        met.add_ensemble(me2)

    def _check_for_clash(self, beads1, beads2):
        count = 0
        for b1 in beads1:
            for b2 in beads2:
                dist = util.distance(b1, b2)
                if dist < settings.CLASH_RADIUS:
                    count += 1
        return count

    def _sample_steps(self, me1, me2):
        mtst = motif_tree_state.MotifTreeStateTree(sterics=0)
        for ms1 in me1.motif_states:
            mtst.add_state(ms1.mts)
            for ms2 in me2.motif_states:
                mtst.add_state(ms2.mts)
                count = self._check_for_clash(mtst.nodes[1].beads,
                                              mtst.nodes[2].beads)

                if count > 2:
                    print ms1.mts.name, ms2.mts.name, count
                mtst.remove_node(mtst.last_node)
            mtst.remove_node(mtst.last_node)

    def test_all_step_ensembles(self):
        bps = ["AU","UA","CG","GC","GU","UG"]
        seen = []
        ensembles = []
        for i,bp1 in enumerate(bps):
            for j,bp2 in enumerate(bps):
                name =  bp1 + "=" + bp2
                if name in seen:
                    continue
                ensembles.append((bp1, bp2,
                                  motif_ensemble.MotifEnsemble(name, 0, 0)))

        for me1 in ensembles:
            for me2 in ensembles:
                if me1[1] != me2[0]:
                    continue
                self._sample_steps(me1[2], me2[2])

    def _get_step_motifs(self, steps):
        motif_names = []
        for i in range(1,len(steps)):
            full_step = steps[i-1] + "=" + steps[i]
            motif_names.append(full_step)
        return motif_names

    def test_helix(self):
        flow_steps = "CG GC AU GC GC AU UA CG UA UA".split()
        motif_names = self._get_step_motifs(flow_steps)
        met = motif_ensemble_tree.MotifEnsembleTree()
        for i, mname in enumerate(motif_names):

            if i == 0:
                met.add_ensemble(motif_ensemble.MotifEnsemble(mname, 0, 0))
            else:
                met.add_ensemble(motif_ensemble.MotifEnsemble(mname, 0, 0))
            #if i == 4:
            #    break
        mtst = met.get_mtst()
        for i, n in enumerate(mtst.nodes):
            if i == 0:
                continue
            print i, n.mts.name
            print n.mts.end_states[1].r
            basic_io.points_to_pdb("beads."+str(i)+".pdb", n.beads)
        mtst.nodes_to_pdbs()
        #mtst.to_pdb("test2.pdb")

        mt = motif_tree.MotifTree(sterics=0)
        for i, n in enumerate(mtst.nodes):
            if i == 0:
                continue
            m = motif.str_to_motif(n.mts.build_string)
            if i == 1:
                mt.add_motif(m, end_index=0, end_flip=0)
            else:
                mt.add_motif(m, end_index=0, end_flip=0)
        mt.write_pdbs("mt_node")

    def test_1(self):
        me = motif_ensemble.MotifEnsemble("GC=GC")
        states = []
        for ms in me.motif_states:
            name_elements = motif_tree_state.parse_db_name(ms.mts.name)
            if name_elements.motif_name == "GC=GC.4":
                states.append(ms.mts)
        for mts in states:
            print mts.name
            for s in mts.end_states:
                if s is not None:
                    print s.r
                    print s.d

    def test_2(self):
        me = motif_ensemble.MotifEnsemble("GC=GC")
        ms = me.get_state("GC=GC.4-0-0-0-0-1-0")
        m = motif.str_to_motif(ms.mts.build_string)
        mt = motif_tree.MotifTree()
        for i in (0, 1):
            for j in (0, 1):
                node = mt.add_motif(m, end_index=i, end_flip=j)
                print i,j,node.available_ends()[0].r()
                node.motif.to_pdb("motif."+str(i)+"."+str(j)+".pdb")
                mt.remove_node(node)

        states = []
        for ms in me.motif_states:
            name_elements = motif_tree_state.parse_db_name(ms.mts.name)
            if name_elements.motif_name == "GC=AU.2":
                states.append(ms.mts)
        mtst = motif_tree_state.MotifTreeStateTree()
        for mts in states:
            mtst.add_state(mts)
            mt = mtst.to_motiftree(sterics=0)
            name_elements = motif_tree_state.parse_db_name(mts.name)
            mt.nodes[1].motif.to_pdb("node."+str(name_elements.start_index)+\
                                     "."+str(name_elements.flip_direction)+".pdb")
            mtst.remove_node(mtst.last_node)

            met.add_ensemble(motif_ensemble.MotifEnsemble(mname, 0, 1))

    def test_3(self):
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        for ms in me.motif_states:
            for s in ms.mts.end_states:
                if s is not None:
                    if s.r[1][0] > 0:
                        print ms.mts.name
                        print s.r

    def test_4(self):
        me = motif_ensemble.MotifEnsemble("GC=GC", 0, 0)
        mtst = motif_tree_state.MotifTreeStateTree(sterics=0)
        mtst.add_state(me.motif_states[0].mts)
        mtst.add_state(me.get_state("GC=GC.100-0-0-0-0-1-0").mts)
        mtst.add_state(me.motif_states[0].mts)
        mtst.nodes_to_pdbs()

    def test_ss_to_ensemble(self):
        sstree = ss_tree.SecondaryStructureTree("(((())))","ACGCGCGU")










def main():
    unittest.main()

if __name__ == '__main__':
    main()
