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

def get_twoway_helix_mts_tree(size=2):
    twoways = motif_tree_state.MotifTreeStateLibrary(motif_type.TWOWAY)
    helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
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

class BuildTaskAstarUnittest(unittest.TestCase):

    def test_creation(self):
        pass

    def test_selector(self):
        pass

    def test_search(self):
        mtss = build_task_astar.MotifTreeStateSearch()

    def test_search_search(self):
        return
        selector =build_task_astar.MotifTreeStateSelector([motif_type.TWOWAY],"all")
        mtss = build_task_astar.MotifTreeStateSearch(max_node_level=2,
                                                     max_solutions=1,
                                                     accept_score=1)
        mtst = get_twoway_mts_tree()
        start = mtst.nodes[0].active_states()[0]
        end = mtst.nodes[-1].active_states()[0]
        solutions = mtss.search(start, end, selector)
        for i in range(len(mtst.nodes)):
            if mtst.nodes[i].mts.name != solutions[0].path[i].mts.name:
                self.fail("did not get correct path back")

    def test_search_search_2(self):
        return
        mtst = get_twoway_helix_mts_tree()
        start = mtst.nodes[0].active_states()[0]
        end = mtst.nodes[-1].active_states()[0]
        mtss = build_task_astar.MotifTreeStateSearch(max_node_level=2,
                                                     max_solutions=1,
                                                     accept_score=1)
        solutions = mtss.search(start, end)
        for i in range(len(mtst.nodes)):
            if mtst.nodes[i].mts.name != solutions[0].path[i].mts.name:
                self.fail("did not get correct path back")

    def test_search_solutions(self):
        mtss =  build_task_astar.MotifTreeStateSearch(max_solutions=1)
        mtst = get_twoway_mts_tree(size=10)
        mtst.to_pdb()
        start = mtst.nodes[0].active_states()[0]
        end = mtst.nodes[-1].active_states()[0]
        solutions = mtss.search(start, end)
        mtst = solutions[0].to_mtst()
        mtst.nodes_to_pdbs()




def main():
    unittest.main()

if __name__ == '__main__':
    main()
