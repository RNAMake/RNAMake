import unittest
import build
import rnamake.motif_state_search

class MotifStateSearchUnittest(unittest.TestCase):

    def test_creation(self):
        mss = rnamake.motif_state_search.MotifStateSearch()

    def test_search(self):
        builder = build.BuildMotifTree()
        mt = builder.build(6)
        mt.write_pdbs("org")
        start = mt.get_node(0).data.ends[0].state()
        end   = mt.get_node(5).data.ends[1].state()
        mss = rnamake.motif_state_search.MotifStateSearch()
        mss.constraint('max_node_level', 5)
        mss.constraint('max_solutions', 1)
        solutions = mss.search(start, end)
        mst = solutions[0].to_mst()
        mst.write_pdbs()

def main():
    unittest.main()

if __name__ == '__main__':
    main()