import sys
import heapq
import motif_tree_state
import motif_type
import base

class PriorityQueue(object):
    def __init__(self):
        self.elements = []

    def empty(self):
        return len(self.elements) == 0

    def put(self, item, priority):
        heapq.heappush(self.elements, (priority, item))

    def get(self):
        return heapq.heappop(self.elements)[1]

class MotifTreeStateSelector(object):
    def __init__(self, motif_types, mode="helix_flank"):
        self.mts_libs = []
        for mtype in motif_types:
            mts_lib = motif_tree_state.MotifTreeStateLibrary(mtype)
            self.mts_libs.append(mts_lib)

    def get_children_mts(self, node):
        pass

class MotifTreeStateSearchScorer(object):
    def __init__(self, target):
        pass

class MotifTreeStateSearch(base.Base):
    def __init__(self, **options):
        self.queue = PriorityQueue()
        self.aligner = motif_tree_state.MotifTreeStateNodeAligner()
        self.scorer = None
        self.lookup = None
        self.steps = 0
        self.solutions = []


