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

    def _test_creation(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)

        mst = motif_state_tree.MotifStateTree()

        for i, n in enumerate(mt):
            if i == 0:
                mst.add_state(rm.manager.get_state(n.data.name))
            else:
                parent_index = 1000
                parent_end_index = -1
                for c in n.connections:
                    if c is None:
                        continue
                    if c.partner(n.index).index < parent_index:
                        parent_index =  c.partner(n.index).index
                        parent_end_index = c.end_index(c.partner(n.index).index)

                if parent_index == 1000:
                    raise ValueError("did not convert motif tree to motif state tree properly")
                mst.add_state(rm.manager.get_state(n.data.name), parent_index, parent_end_index)

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

    def test_res_ids(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        mst = motif_state_tree.MotifStateTree(mt)

        for r in mt.residues():
            r_mst = mst.get_residue(r.uuid)
            if r_mst is None:
                raise ValueError("cannot find residue in mst")


    def test_align(self):
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

    def _test_replace_state(self):
        builder = build.BuildMotifTree()
        mt = builder.build(10)
        mst = motif_state_tree.MotifStateTree(mt)
        rstate = rm.manager.get_state("HELIX.IDEAL.2")
        mst.replace_state(4, rstate)

        mt2 = motif_tree.MotifTree()
        for n in mst:
            m = rm.manager.get_motif(n.data.ref_state.name)
            mt2.add_motif(m)

        i = len(mst)-1
        d1 = mst.get_node(i).data.cur_state.end_states[1].d
        d2 = mt2.get_node(i).data.ends[1].d()

        if util.distance(d1, d2) > 0.1:
            self.fail("did not properly replace state")

    def _test_specific(self):
        motifs_str = """HELIX.IDEAL.16
TWOWAY.3R1C.5
HELIX.IDEAL.7
TWOWAY.3DIL.0
HELIX.IDEAL.12
TWOWAY.3LOA.0
HELIX.IDEAL.11
TWOWAY.2VQE.3
HELIX.IDEAL.14
TWOWAY.4K27.0"""
        motifs = motifs_str.split("\n")
        builder = build.BuildMotifTree()
        mt = builder.build_specific(motifs)
        mt.write_pdbs()
        mst = motif_state_tree.MotifStateTree(mt)
        i = 0
        for n in mst:
            basic_io.points_to_pdb("beads."+str(i)+".pdb",n.data.cur_state.beads)
            i += 1








def main():
    unittest.main()

if __name__ == '__main__':
    main()