import sys
import heapq
import motif_tree_state
import motif_type
import motif
import base
import option
import util
import math
import transform
import settings
import logger
import motif_tree_state_selector
import numpy as np


clogger =  logger.get_logger("MotifTreeStateSearch:search")


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


class MotifTreeStateSearchNode(motif_tree_state.MotifTreeStateNode):
    def __init__(self, mts, parent, lib_type):
        if parent is None:
            self.level = 0
            self.node_counts = []
        else:
            self.level = parent.level+1
        self.lib_type = lib_type
        self.mts, self.parent = mts, parent
        self.beads, self.score, self.size, self.ss_score = mts.beads, 1000, 0, 0
        self.states = [None for s in self.mts.end_states]
        for i, s in enumerate(self.mts.end_states):
            if s is None:
                continue
            self.states[i] = s.copy()

    def copy(self):
        c = MotifTreeStateSearchNode(self.mts, self.parent, self.lib_type)
        c.beads = np.copy(self.beads)
        c.score, c.size, c.ss_score = self.score, self.size, c.ss_score
        c.states = [None for s in self.states]
        c.level = self.level
        for i, s in enumerate(self.states):
            if s is None:
                continue
            c.states[i] = s.copy()
        return c

    def lib_type_usage(self, lib_type):
        return self.node_counts[lib_type]

    def update_node_count(self):
        self.node_counts = [x for x in self.parent.node_counts]
        self.node_counts[self.lib_type] += 1


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
        self.node_selector = motif_tree_state_selector.default_selector()
        self.lookup = None
        self.using_lookup = 0
        self.steps = 0
        self.solutions = []
        self.setup_options_and_constraints()
        self._set_option_or_constraint(options)

    def search(self, start, end, node_selector=None, lookup=None, **options):
        if node_selector is not None:
            self.node_selector = node_selector
        if lookup is not None:
            self.lookup = lookup
            self.using_lookup = 1
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
        verbose        = self.option('verbose')
        frequency      = self.option('frequency')
        while not self.queue.empty():
            current = self.queue.get()
            if current.level > max_node_level:
                continue
            self.steps += 1
            if verbose and self.steps % frequency == 0:
                self._print_status()

            score = self.scorer.accept_score(current)
            if score < accept_score:
                fail=0
                if not self.node_selector.is_valid_solution(current):
                    fail = 1
                if sterics and not fail:
                    if current.steric_clash():
                        fail = 1
                if not fail:
                    solution = MotifTreeStateSearchSolution(current, score)
                    self.solutions.append(solution)
                    if len(self.solutions) == self.constraint('max_solutions'):
                        return self.solutions

            if current.level == max_node_level:
                continue

            children, types = self.node_selector.get_children_mts(current)
            test_node.parent = current
            parent_ends = current.active_states()
            test_node.level = current.level+1

            for i, c in enumerate(children):
                for end in parent_ends:
                    test_node.replace_mts(c)
                    self.aligner.transform_state(end, current, test_node)
                    score = self.scorer.score(test_node)
                    if score > current.score:
                        continue
                    if sterics:
                        self.aligner.transform_beads(test_node)
                        if test_node.steric_clash(2):
                            continue
                        if self.using_lookup:
                            if self.lookup.clash(test_node.beads):
                                continue


                    child = test_node.copy()
                    child.lib_type = types[i]
                    child.update_node_count()
                    self.queue.put(child, score)

        return self.solutions

    def setup_options_and_constraints(self):
        options =     { 'sterics'        :  1,
                        'verbose'        :  0,
                        'frequency'      : 10 }

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
                continue
            except:
                pass

            try:
                self.constraints.set(k, v)
            except:
                raise ValueError("cannot set "+k+" with value "+str(v)+ \
                                 "it is neither an option nor a constraint")

    def _get_start_node(self, start):
        mts = motif_tree_state.MotifTreeState("start", 1, 0, 0, [],
                                              [start], 0, "")
        start_node = MotifTreeStateSearchNode(mts, None, -1)
        start_node.node_counts = [0 for x in self.node_selector.nodes]

        return start_node

    def _print_status(self):
        clogger.info("steps: "+ str(self.steps))
    def reset(self):
        self.queue = PriorityQueue()
        self.solutions = []
        self.steps = 0


class MotifTreeStateSearchSolution(object):
    def __init__(self, node, score):
        self.path = self._get_path(node)

    def _get_path(self, node):
        path = []
        while node != None:
            path.append(node)
            node = node.parent
        return path[::-1]

    def to_mtst(self):
        start_mts = self.path[0].mts
        ref_motif = motif.ref_motif()
        r = util.unitarize(start_mts.end_states[0].r.T.dot(ref_motif.ends[0].r()))
        trans = -ref_motif.ends[0].d()
        t = transform.Transform(r, trans)
        ref_motif.transform(t)
        bp_pos_diff = start_mts.end_states[0].d - ref_motif.ends[0].d()
        ref_motif.move(bp_pos_diff)
        start_mts.build_string = ref_motif.to_str()
        start_mts.name = "mtss_start"
        mtst = motif_tree_state.MotifTreeStateTree(start_mts, sterics=0)
        for i, n in enumerate(self.path):
            if i == 0:
                continue
            success = 0
            for s in mtst.last_node.active_states():
                node = mtst.add_state(n.mts, parent_end=s)
                if node is None:
                    continue
                same=1
                for j in range(len(node.states)):
                    if n.states[j] is None and node.states[j] is None:
                        continue
                    if n.states[j] is None and node.states[j] is not None:
                        same=0
                        break
                    if n.states[j] is not None and node.states[j] is None:
                        same=0
                        break
                    diff = n.states[j].diff(node.states[j])
                    if diff > 0.1:
                        same=0
                        break
                if same:
                    success=1
                    break
                else:
                    mtst.remove_node(node)
            if not success:
                print "fail"
                exit()

        return mtst

    def nodes_to_pdbs(self):
        self.to_mtst().nodes_to_pdbs()

    def to_pdb(self, fname="mtss_solution.pdb"):
        self.to_mtst().to_pdb(fname)

    def to_pose(self):
        return self.to_mtst().to_pose()

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

