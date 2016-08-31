import motif
import motif_state_tree
import base
import option
import priority_queue
import motif_state_selector
import motif_state_search_scorer
import copy
import numpy as np

class MotifStateSearch(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()
        self.queue = priority_queue.PriorityQueue()
        self.selector = motif_state_selector.default_selector()
        self.scorer = motif_state_search_scorer.MTSS_Astar()
        self.solutions = []
        self.lookup = None
        self.test_node = None
        self.no_more_solutions = 0

    def setup_options_and_constraints(self):
        options =     { 'sterics'        :  1,
                        'verbose'        :  0,
                        'frequency'      :  10,
                        'max_node_level'  : 20,
                        'max_steps'       : 100000000,
                        'max_solutions'   : 10,
                        'max_size'        : 100000000,
                        'min_size'        : 0,
                        'accept_score'    : 10}

        self.options = option.Options(options)

    def _start_node(self, start_bp):
        ms = motif.MotifState('start', ['start', 'start'], ['', ''],
                              [start_bp, start_bp], [], 0, 0, 0)
        n = MotifStateSearchNode(ms, None, -1, -1)
        n.node_type_usages = [0 for i in range(len(self.selector.graph))]
        return n

    def _get_local_variables(self):
        return self.option('accept_score'), self.option('max_node_level'), \
               self.option('sterics'), self.option('max_size'), self.option('min_size')

    def setup(self, start, end):
        start_n = self._start_node(start)
        start_n.score = 10000000
        test_node = start_n.copy()
        self.queue.put(start_n, 10000000)
        self.scorer.set_target(end)
        self.test_node = start_n.copy()

    def finished(self):
        if len(self.solutions) >= self.option("max_solutions") or \
            self.no_more_solutions == 1:
            return 1
        else:
            return 0

    def next(self):
        sol = self._search()
        if sol == None:
            self.no_more_solutions = 1
            return None
        self.solutions.append(sol)
        return sol

    def all(self):
        while not self.queue.empty():
            sol = self._search()
            self.solutions.append(sol)
            if len(self.solutions) >= self.option('max_solutions'):
                break
        return self.solutions

    def _search(self):
        accept_score, max_node_level, sterics, max_size, min_size = self._get_local_variables()
        best = 1000000000
        while not self.queue.empty():
            current = self.queue.get()
            score = self.scorer.accept_score(current)

            #print current.score
            if current.score < best:
                best = current.score
                #print current.score, current.level


            if score < accept_score:
                if not self.selector.is_valid_solution(current):
                    continue
                if current.size < min_size:
                    continue
                return MotifStateSearchSolution(current, score)

            if current.level+1 > max_node_level:
                continue

            motif_states, types = self.selector.get_children_ms(current)
            avail_ends = []
            for i, end in enumerate(current.cur_state.end_states):
                if i == 0:
                    continue
                avail_ends.append(end)
            self.test_node.parent = current
            self.test_node.update()
            for end in avail_ends:
                parent_end_index = current.cur_state.end_states.index(end)
                self.test_node.parent_end_index = parent_end_index
                for i, ms in enumerate(motif_states):
                    self.test_node.update_ms(ms)
                    self.test_node.ntype = types[i]
                    motif.get_aligned_motif_state(end,
                                                  self.test_node.cur_state,
                                                  self.test_node.ref_state)

                    score = self.scorer.score(self.test_node)
                    #score += self.selector.score(self.test_node)*self.test_node.level*10
                    #print score, current.score, len(current.cur_state.beads)
                    if score > current.score:
                        continue

                    if sterics:
                        if self.lookup is not None:
                            if self.lookup.clash(self.test_node.cur_state.beads):
                                continue
                    #print "test", len(self.test_node.cur_state.beads)
                    child = self.test_node.copy()
                    #print "made it"
                    child.score = score
                    child.update()
                    #print "child", len(child.cur_state.beads)
                    #exit()
                    if child.size > max_size:
                        continue
                    child.ntype = types[i]
                    self.queue.put(child, score)


        return None


class MotifStateSearchNode(object):
    def __init__(self, ref_state, parent, parent_end_index, ntype):
        self.parent = parent
        self.parent_end_index = parent_end_index
        self.ss_score = ref_state.score
        self.size =   ref_state.size
        self.node_type_usages = []
        if parent is None:
            self.level = 1
        else:
            self.level = self.parent.level + 1
            self.ss_score += self.parent.ss_score
        self.ref_state, self.parent, self.ntype = ref_state, parent, ntype
        self.cur_state = self.ref_state.copy()
        self.score = 1000

    def copy(self):
        new_n = MotifStateSearchNode(self.ref_state, self.parent,
                                     self.parent_end_index, self.ntype)
        #print "copied", len(self.cur_state.copy().beads)
        new_n.cur_state = self.cur_state.copy()
        new_n.score = self.score
        new_n.ss_score = self.ss_score
        new_n.level = self.level
        new_n.node_type_usages = self.node_type_usages[::]
        return new_n

    def node_type_usage(self, i):
        if i == -1:
            return 0
        return self.node_type_usages[i]

    def update(self):
        if self.parent is None:
            return
        self.level = self.parent.level + 1
        self.ss_score = self.ref_state.score + self.parent.ss_score
        self.size = self.ref_state.size + self.parent.size
        if self.ntype == -1:
            return
        self.node_type_usages = self.parent.node_type_usages[::]
        self.node_type_usages[self.ntype] +=1

    def update_ms(self, ms):
        self.ref_state = ms
        self.cur_state.end_states = [e.copy() for e in ms.end_states]


class MotifStateSearchSolution(object):
    def __init__(self, node, score):
        self.path = self._get_path(node)
        self.score = score

    def _get_path(self, node):
        path = []
        while node != None:
            path.append(node)
            node = node.parent
        return path[::-1][1:]

    def to_mst(self):
        mst = motif_state_tree.MotifStateTree()
        mst.option('sterics', 0)
        for i, n in enumerate(self.path):
            n.cur_state.name = n.ref_state.name
            n.cur_state.end_ids = n.ref_state.end_ids
            n.cur_state.end_names = n.ref_state.end_names
            if i == 0:
                mst.add_state(n.cur_state)
            else:
                j = mst.add_state(n.cur_state, -1, n.parent_end_index)
                if j == -1:
                    raise ValueError("something went horribly wrong, cannot build solution")
        return mst

    def to_motif_tree(self):
        return self.to_mst().to_motif_tree()

