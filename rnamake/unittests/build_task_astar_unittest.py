import unittest
import rnamake.build_task_astar as build_task_astar
import rnamake.motif_type as motif_type
import rnamake.motif_tree_state as motif_tree_state
import rnamake.resource_manager as resource_manger
import rnamake.steric_lookup as steric_lookup
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


def get_nway_helix_mts_tree(size=2):
    nways = motif_tree_state.MotifTreeStateLibrary(motif_type.NWAY)
    helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
    me_libs = [helixs, nways]
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

    def test_search_search_3(self):
        """nways = motif_tree_state.MotifTreeStateLibrary(motif_type.NWAY)
        helixs = motif_tree_state.MotifTreeStateLibrary(motif_type.HELIX)
        me_libs = [helixs, nways]
        mtst = motif_tree_state.MotifTreeStateTree()
        mtst.add_state(helixs.get_state('HELIX.LE.4-0-0-1-4-0-1'))
        mtst.add_state(nways.get_state('NWAY.2VQE.5-0-0-1-0-2-1'))
        mtst.nodes_to_pdbs()"""
        mtst = get_nway_helix_mts_tree(size=4)
        mtst.to_pdb()
        start = mtst.nodes[0].active_states()[0]
        end = mtst.nodes[-1].active_states()[0]
        mtss = build_task_astar.MotifTreeStateSearch(max_solutions=1,
                                                     accept_score=14,
                                                     verbose=1,
                                                     frequency=10)

        ns = build_task_astar.MotifTreeStateSelector([motif_type.NWAY])
        solutions = mtss.search(start, end, ns)
        solutions[0].nodes_to_pdbs()
        for n in solutions[0].path:
            print n.mts.name



    def test_search_solutions(self):
        return
        mtss =  build_task_astar.MotifTreeStateSearch(max_solutions=1)
        mtst = get_twoway_mts_tree(size=10)
        mtst.to_pdb()
        start = mtst.nodes[0].active_states()[0]
        end = mtst.nodes[-1].active_states()[0]
        solutions = mtss.search(start, end)
        mtst = solutions[0].to_mtst()
        mtst.nodes_to_pdbs()

    def test_lookup(self):
        return
        mtss =  build_task_astar.MotifTreeStateSearch(max_solutions=1,
                                                      accept_score=10)
        mtst = get_twoway_mts_tree(size=10)
        p = mtst.to_pose()
        p.to_pdb()
        start = mtst.nodes[0].active_states()[0]
        end = mtst.nodes[-1].active_states()[0]
        sl = steric_lookup.StericLookup(p)
        solutions = mtss.search(start, end, lookup=sl)
        mtst = solutions[0].to_mtst()
        mtst.nodes_to_pdbs()

    def test_search_solutions_non_ideal(self):
        return
        mtss =  build_task_astar.MotifTreeStateSearch(max_solutions=1,
                                                      accept_score=10)
        mtst = get_twoway_mts_tree(size=4)
        mtst.to_pdb()
        s = None
        for state in mtst.nodes[2].states:
            if state is not None:
                s = state
                break
        end = mtst.nodes[-1].active_states()[0]
        solutions = mtss.search(s, end)
        mtst = solutions[0].to_mtst()
        mtst.nodes_to_pdbs()






def main():
    unittest.main()

if __name__ == '__main__':
    main()
