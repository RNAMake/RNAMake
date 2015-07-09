import unittest
import build
import rnamake.motif_state_search
import rnamake.util as util
import rnamake.motif_state_tree as motif_state_tree
import rnamake.settings as settings
import rnamake.resource_manager as rm
import rnamake.steric_lookup as steric_lookup
import rnamake.motif_state_search_scorer as motif_state_search_scorer

class MotifStateSearchUnittest(unittest.TestCase):

    def test_creation(self):
        mss = rnamake.motif_state_search.MotifStateSearch()

    def test_search(self):
        builder = build.BuildMotifTree()
        mt = builder.build(3)
        start = mt.get_node(0).data.ends[0].state()
        end   = mt.last_node().data.ends[1].state()
        mss = rnamake.motif_state_search.MotifStateSearch()
        mss.constraint('max_node_level', 3)
        mss.constraint('max_solutions', 1)
        solutions = mss.search(start, end)
        mst = solutions[0].to_mst()
        new_end = mst.last_node().data.cur_state.end_states[1]

        dist = util.distance(new_end.d, end.d)
        if dist > mss.constraint('accept_score'):
            self.fail("did not find a suitable solution")

    def test_search_2(self):
        path = settings.UNITTEST_PATH + "/resources/motifs/tetraloop_receptor_min"
        rm.manager.add_motif(path)
        mst = motif_state_tree.MotifStateTree(sterics=0)
        ms = rm.manager.get_state("tetraloop_receptor_min", "A228-A246")
        m  = rm.manager.get_motif("tetraloop_receptor_min", "A228-A246")
        mst.add_state(ms)
        start = mst.get_node(0).data.get_end_state("A221-A252")
        end   = mst.get_node(0).data.get_end_state("A146-A157")
        mss = rnamake.motif_state_search.MotifStateSearch()
        mss.constraint('max_node_level', 6)
        mss.constraint('max_solutions', 1)
        mss.constraint('accept_score', 10)
        mss.scorer = motif_state_search_scorer.MTSS_Astar()

        sl = steric_lookup.StericLookup()
        beads = m.get_beads(m.ends)
        sl.add_beads(beads)
        mss.lookup = sl

        solutions = mss.search(start, end)
        mst_sol = solutions[0].to_mst()
        mst.add_mst(mst_sol, parent_end_name="A221-A252")
        mst.add_connection(0, mst.last_node().index, "A146-A157")








def main():
    unittest.main()

if __name__ == '__main__':
    main()