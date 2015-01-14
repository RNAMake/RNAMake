import unittest
import rnamake.build_task_astar
import rnamake.motif_type as motif_type
import rnamake.motif_tree_state as motif_tree_state
import rnamake.resource_manager as resource_manger

def get_twoway_mts_tree(self, size=2):
    pass

class BuildTaskAstarUnittest(unittest.TestCase):

    def test_creation(self):
        pass

    def test_selector(self):
        selector = rnamake.build_task_astar.MotifTreeStateSelector([motif_type.TWOWAY])

    def test_search(self):
        mtss = rnamake.build_task_astar.MotifTreeStateSearch()

    def test_search_search(self):
        pass


def main():
    unittest.main()

if __name__ == '__main__':
    main()
