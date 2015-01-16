import unittest
import rnamake.build_task_astar as build_task_astar
import rnamake.motif_type as motif_type
import rnamake.motif_tree_state as motif_tree_state
import rnamake.resource_manager as resource_manger
import random

def get_twoway_mts_tree(size=2):
    mts_lib = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
    mtst = motif_tree_state.MotifTreeStateTree()
    while len(mtst.nodes) < size+1:
        mts = random.choice(mts_lib.motif_tree_states)
        mtst.add_state(mts)

    return mtst


class BuildTaskAstarUnittest(unittest.TestCase):

    def test_creation(self):
        pass

    def test_selector(self):
        pass

    def test_search(self):
        mtss = build_task_astar.MotifTreeStateSearch()

    def test_search_search(self):
        return
        mtss = build_task_astar.MotifTreeStateSearch(max_node_level=2,
                                                     max_solutions=1,
                                                     accept_score=1)
        mtst = get_twoway_mts_tree()
        start = mtst.nodes[0].active_states()[0]
        end = mtst.nodes[-1].active_states()[0]
        solutions = mtss.search(start, end)
        for i in range(len(mtst.nodes)):
            if mtst.nodes[i].mts.name != solutions[0].path[i].mts.name:
                self.fail("did not get correct path back")


    def test_search_solutions(self):
        mtss =  build_task_astar.MotifTreeStateSearch(max_solutions=1)
        mtst = get_twoway_mts_tree(size=10)
        start = mtst.nodes[0].active_states()[0]
        end = mtst.nodes[-1].active_states()[0]
        solutions




def main():
    unittest.main()

if __name__ == '__main__':
    main()
