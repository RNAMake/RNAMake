import sys
import heapq
import motif_tree_state
import motif_type

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
    def __init__(self, motif_types):
        pass
