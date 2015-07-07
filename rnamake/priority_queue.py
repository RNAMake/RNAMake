import heapq

class PriorityQueue(object):
    """
    Helper class for deciding which nodes to visit first during
    MotifTreeStateSearch. Organizes nodes by score using heapq modulue
    """

    def __init__(self):
        self.elements = []

    def empty(self):
        """
        return true if there are no elements false if there is
        """
        return len(self.elements) == 0

    def put(self, item, priority):
        """
        adds new element with a given score to the queue
        """
        heapq.heappush(self.elements, (priority, item))

    def get(self):
        """
        gets next best element
        """
        return heapq.heappop(self.elements)[1]

