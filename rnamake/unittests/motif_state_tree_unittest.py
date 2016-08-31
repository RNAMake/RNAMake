import unittest
import rnamake.motif as motif
import rnamake.motif_tree as motif_tree
import rnamake.motif_factory as motif_factory
import rnamake.motif_state_tree as motif_state_tree
import rnamake.resource_manager as rm
import rnamake.util as util
import rnamake.settings as settings
import rnamake.basic_io as basic_io
import rnamake.eternabot.sequence_designer as sequence_designer
import build

class MotifStateTreeUnittest(unittest.TestCase):


    def test_creation_from_mt(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(mt)
        if len(mst) != 10:
            self.fail("did not build mst properly")

        #for n in mt:
        #    print n.data.name

        #print mt.last_node().data.ends[0].d()
        #print mst.last_node().data.cur_state.end_states[0].d

    #TODO fix
    def _test_align(self):
        path = settings.UNITTEST_PATH + "/resources/motifs/tetraloop_receptor_min"
        rm.manager.add_motif(path)
        m  = rm.manager.get_motif(name="tetraloop_receptor_min", end_name="A228-A246")
        bp_state = m.ends[1].state()
        test_state = rm.manager.ms_libs["ideal_helices"].get(name='HELIX.IDEAL.3')
        d1 = bp_state.d
        #rint d1
        motif.align_motif_state(bp_state, test_state)
        d2 = test_state.end_states[0].d
        if util.distance(d1, d2) > 0.5:
            self.fail("did not align motif state properly")

    def _test_change_sequence(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        mst = motif_state_tree.MotifStateTree(mt)
        ss = mst.designable_secondary_structure()
        designer = sequence_designer.SequenceDesigner()
        results = designer.design(ss.dot_bracket(), ss.sequence())
        ss.replace_sequence(results[0].sequence)
        connectivity = ss.motif_topology_from_end(ss.ends[0])

    def _test_topology_to_str(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        mst = motif_state_tree.MotifStateTree(mt)
        s = mst.topology_to_str()
        mst2 = motif_state_tree.str_to_motif_state_tree(s)
        if len(mst2) != len(mt):
            self.fail("did not reconstituate motif_state_tree from topology")

    def test_remove_node(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mst = motif_state_tree.MotifStateTree(mt)
        mst.remove_node(mst.last_node().index)
        if len(mst) != 2:
            self.fail("did not properly remove node")

    def test_remove_node_level(self):
        pass

    def test_replace_state(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(mt)
        rstate = rm.manager.get_state(name="HELIX.IDEAL.2")
        mst.replace_state(2, rstate)

        mt2 = motif_tree.MotifTree()
        for n in mst:
            m = rm.manager.get_motif(name=n.data.cur_state.name,
                                     end_name=n.data.cur_state.end_names[0])
            mt2.add_motif(m)

        i = len(mst)-1
        d1 = mst.get_node(i).data.cur_state.end_states[1].d
        d2 = mt2.get_node(i).data.ends[1].d()
        if util.distance(d1, d2) > 0.1:
            self.fail("did not properly replace state")








def main():
    unittest.main()

if __name__ == '__main__':
    main()