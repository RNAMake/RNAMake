import unittest
import build
import rnamake.motif_state_search
import rnamake.util as util
import rnamake.motif_state_tree as motif_state_tree
import rnamake.settings as settings
import rnamake.resource_manager as rm
import rnamake.steric_lookup as steric_lookup
import rnamake.motif_state_search_scorer as motif_state_search_scorer
import rnamake.segmenter as segmenter
import rnamake.pose_factory as pf
from rnamake import motif_tree

class MotifStateSearchUnittest(unittest.TestCase):

    def test_creation(self):
        mss = rnamake.motif_state_search.MotifStateSearch()

    def test_search(self):
        builder = build.BuildMotifTree()
        mt = builder.build(2)
        start = mt.get_node(0).data.ends[0].state()
        end   = mt.last_node().data.ends[1].state()
        mss = rnamake.motif_state_search.MotifStateSearch()
        mss.option('max_node_level', 3)
        mss.option('accept_score', 2)
        mss.setup(start, end)
        s = mss.next()
        if s is None:
            print mt
            raise ValueError("could not find a suitable solution")
        mst = s.to_mst()
        new_end = mst.last_node().data.cur_state.end_states[1]

        dist = util.distance(new_end.d, end.d)
        if dist > mss.option('accept_score'):
            self.fail("did not find a suitable solution")

    def test_search_2(self):
        path = settings.UNITTEST_PATH + "/resources/motifs/tetraloop_receptor_min"
        rm.manager.add_motif(path)
        mst = motif_state_tree.MotifStateTree(sterics=0)
        ms = rm.manager.get_state(name="tetraloop_receptor_min", end_name="A228-A246")
        m  = rm.manager.get_motif(name="tetraloop_receptor_min", end_name="A228-A246")
        mst.add_state(ms)
        start = mst.get_node(0).data.get_end_state("A221-A252")
        end   = mst.get_node(0).data.get_end_state("A146-A157")
        mss = rnamake.motif_state_search.MotifStateSearch()
        mss.option('max_node_level', 6)
        mss.option('max_solutions', 1)
        mss.scorer = motif_state_search_scorer.MTSS_Astar()

        sl = steric_lookup.StericLookup()
        beads = m.get_beads(m.ends)
        sl.add_beads(beads)
        mss.lookup = sl

        mss.setup(start, end)
        s = mss.next()
        mst_sol = s.to_mst()
        mst.add_mst(mst_sol, parent_end_name="A221-A252")
        mst.add_connection(0, mst.last_node().index, "A146-A157")
        #p = mst.to_pose()

    def test_redesign(self):
        s = rnamake.segmenter.Segmenter()
        path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
        p = pf.factory.pose_from_file(path)
        end1 = p.get_basepair(name='A111-A209')[0]
        end2 = p.get_basepair(name='A118-A203')[0]
        segments = s.apply(p, [end1, end2])
        pd = segments.remaining
        start = end1.state()
        end   = end2.state()

        mss = rnamake.motif_state_search.MotifStateSearch()
        mss.option('max_node_level', 5)
        mss.option('max_solutions', 1)
        mss.option('accept_score', 5)
        mss.scorer = motif_state_search_scorer.MTSS_Astar()

        sl = steric_lookup.StericLookup()
        beads = pd.get_beads(pd.ends)
        sl.add_beads(beads)
        mss.lookup = sl

        mss.setup(start, end)
        s = mss.next()
        mst_sol = s.to_mst()







def main():
    unittest.main()

if __name__ == '__main__':
    main()