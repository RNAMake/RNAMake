import sys
import heapq
import motif_tree_state
import motif_type
import base
import option
import util
import math

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
        return self.mts_libs[0].motif_tree_states


class MotifTreeStateSearchScorer(object):
    def __init__(self, target):
        self.target = target
        self.target_flip = target.copy()
        self.target_flip.flip()

    def score(self, node):
        raise ValueError("cannot call base MotifTreeStateSearchScorer")

    def accept_score(self, node):
        current = node.active_states()[0]
        score = util.distance(current.d, self.target.d)

        r_diff       = util.matrix_distance(current.r, self.target.r)
        r_diff_flip  = util.matrix_distance(current.r, self.target_flip.r)

        if r_diff > r_diff_flip:
            r_diff = r_diff_flip

        score += 2*r_diff
        return score


class MTSS_GreedyBestFirstSearch(MotifTreeStateSearchScorer):
    def __init__(self, target):
        super(self.__class__, self).__init__(target)

    def score(self, node):
        return new_score_function(node.active_states()[0], self.target,
                                  self.target_flip)


class MotifTreeStateSearch(base.Base):
    def __init__(self, **options):
        self.queue = PriorityQueue()
        self.aligner = motif_tree_state.MotifTreeStateNodeAligner()
        self.scorer = None
        self.node_selector = MotifTreeStateSelector([motif_type.TWOWAY])
        self.lookup = None
        self.steps = 0
        self.solutions = []
        self.setup_options_and_constraints()
        self._set_option_or_constraint(options)

    def search(self, start, end, node_selector=None, **options):
        self.scorer = MTSS_GreedyBestFirstSearch(end)
        start_node = self._get_start_node(start)
        test_node = start_node.copy()
        self.queue.put(start_node, 10000)
        children = []
        current = None
        score = 0
        accept_score   = self.constraint('accept_score')
        max_node_level = self.constraint('max_node_level')
        sterics        = self.option('sterics')
        while not self.queue.empty():
            current = self.queue.get()
            if current.level > max_node_level:
                continue


            score = self.scorer.accept_score(current)
            if score < accept_score:
                solution = MotifTreeStateSearchSolution(current, score)
                self.solutions.append(solution)
                if len(self.solutions) == self.constraint('max_solutions'):
                    return self.solutions

            if current.level == max_node_level:
                continue

            children = self.node_selector.get_children_mts(current)
            test_node.parent = current
            parent_ends = current.active_states()
            test_node.level = current.level+1

            for c in children:
                test_node.replace_mts(c)
                self.aligner.transform_state(parent_ends[0], current, test_node)
                score = self.scorer.score(test_node)
                if score > current.score:
                    continue
                if sterics:
                   self.aligner.transform_beads(test_node)

                child = test_node.copy()
                self.queue.put(child, score)

        return self.solutions

    def setup_options_and_constraints(self):
        options =     { 'sterics'        :  1,
                        'vebose'         :  0 }

        constraints = { 'max_node_level'  : 999,
                        'max_steps'       : 100000000,
                        'max_solutions'   : 10,
                        'accept_score'    : 10}

        self.options = option.Options(options)
        self.constraints = option.Options(constraints)

    def _set_option_or_constraint(self, options):
        for k, v in options.iteritems():
            try:
                self.options.set(k, v)
            except:
                pass

            try:
                self.constraints.set(k, v)
            except:
                raise ValueError("cannot set "+k+" with value "+v+" it is " \
                                 "it is neither an option nor a constraint")

    def _get_start_node(self, start):
        mts = motif_tree_state.MotifTreeState("start", 1, 0, 0, [],
                                              [start], 0, "")
        start_node = motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])
        return start_node


class MotifTreeStateSearchSolution(object):
    def __init__(self, node, score):
        self.path = self._get_path(node)

    def _get_path(self, node):
        path = []
        while node != None:
            path.append(node)
            node = node.parent
        return path[::-1]



def new_score_function(current, end, endflip):
    d_diff = util.distance(current.d,end.d)

    if d_diff > 25:
        return d_diff

    r_diff       = util.matrix_distance(current.r, end.r)
    r_diff_flip  = util.matrix_distance(current.r, endflip.r)

    if r_diff > r_diff_flip:
        r_diff = r_diff_flip

    if d_diff < 0.0001:
        d_diff = 0.00001
    scale = (math.log(150/d_diff) - 1)
    if scale > 2:
        scale = 2
    if scale < 0:
        scale = 0

    return d_diff + scale*r_diff

