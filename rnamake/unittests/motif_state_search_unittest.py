import unittest
import build

from rnamake import motif_state_search, resource_manager, motif_tree, util

class MotifStateSearchUnittest(unittest.TestCase):

    def setUp(self):
        self.rm = resource_manager.ResourceManager()


    def _test_creation(self):
        self.mss = motif_state_search.MotifStateSearch(self.rm)

    def _test_search(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name='HELIX.IDEAL.10')
        mt.add_motif(m_name='TWOWAY.2VQE.19', m_end_name='A1008-A1021')

        start = mt.get_node(0).data
        end   = mt.last_node().data
        mss = motif_state_search.MotifStateSearch(self.rm)
        mss.option('max_node_level', 3)
        mss.option('accept_score', 0.5)
        mss.setup(start, 0,  end, 1)
        s = mss.next()
        self.failUnless(s is not None)
        mst = s.to_mst(self.rm)
        new_end = mst.last_node().data.get_end(1)

        dist = util.distance(new_end.d, mt.last_node().data.get_end(1).d)
        if dist > mss.option('accept_score'):
            self.fail("did not find a suitable solution")

    def test_search_longer(self):
        mt = motif_tree.MotifTree(self.rm)
        mt.add_motif(m_name='HELIX.IDEAL.10')
        mt.add_motif(m_name='TWOWAY.2VQE.19', m_end_name='A1008-A1021')
        mt.add_motif(m_name='HELIX.IDEAL.10')
        mt.add_motif(m_name='TWOWAY.2VQE.19', m_end_name='A1008-A1021')
        mt.add_motif(m_name='HELIX.IDEAL.10')
        mt.add_motif(m_name='TWOWAY.2VQE.19', m_end_name='A1008-A1021')
        mt.add_motif(m_name='HELIX.IDEAL.10')
        mt.add_motif(m_name='TWOWAY.2VQE.19', m_end_name='A1008-A1021')
        mt.to_pdb("test.pdb", renumber=1)

        start = mt.get_node(0).data
        end = mt.last_node().data
        mss = motif_state_search.MotifStateSearch(self.rm)
        mss.option('max_node_level', 10)
        mss.setup(start, 0, end, 1)
        s = mss.next()

        mt2 = s.to_motif_tree(self.rm)
        mt2.to_pdb("search.pdb", renumber=1)

    #TODO fix
    """def _test_search_2(self):
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

    def _test_redesign(self):
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
        mss.scorer = motif_state_search_scorer.MTSS_Astar()o

        sl = steric_lookup.StericLookup()
        beads = pd.get_beads(pd.ends)
        sl.add_beads(beads)
        mss.lookup = sl

        mss.setup(start, end)
        s = mss.next()
        mst_sol = s.to_mst()

        """





def main():
    unittest.main()

if __name__ == '__main__':
    main()
