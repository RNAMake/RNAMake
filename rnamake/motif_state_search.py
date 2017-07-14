import motif
import motif_state_tree
import base
import option
import priority_queue
import motif_state
import motif_state_selector
import motif_state_search_scorer

from collections import defaultdict

class MotifStateSearch(base.Base):
    __slots__ = [
        '_rm',
        '_queue',
        '_selector',
        '_scorer',
        '_solutions',
        '_test_node',
        '_no_more_solutions',
        'options'
    ]

    def __init__(self, rm):
        self._rm = rm
        self._queue = priority_queue.PriorityQueue()
        self._selector = motif_state_selector.default_selector(self._rm)
        self._scorer = motif_state_search_scorer.MTSS_Astar()
        self._solutions = []

        self._no_more_solutions = 0
        self.setup_options_and_constraints()
        #self.lookup = None
        #self.test_node = None

    def setup_options_and_constraints(self):
        options =     { 'sterics'        :  1,
                        'verbose'        :  0,
                        'frequency'      :  10,
                        'max_node_level'  : 20,
                        'max_steps'       : 100000000,
                        'max_solutions'   : 1,
                        'max_size'        : 100000000,
                        'min_size'        : 0,
                        'accept_score'    : 10.0}

        self.options = option.Options(options)

    def _get_local_variables(self):
        return self.option('accept_score'), self.option('max_node_level'), \
               self.option('sterics'), self.option('max_size'), self.option('min_size')

    def setup(self, start, start_ei, end, end_ei):
        start_n = MotifStateSearchNode(start, None, -1, -1, 10000000)
        self._test_node = MotifStateSearchNode(start, None, -1, -1, 1000000)
        self._scorer.set_target(end.get_end(end_ei))
        self._first_round(start_n, start_ei)

    def finished(self):
        if len(self._solutions) >= self.option("max_solutions") or \
            self.no_more_solutions == 1:
            return 1
        else:
            return 0

    def next(self):
        sol = self._search()
        if sol == None:
            self.no_more_solutions = 1
            return None
        self._solutions.append(sol)
        return sol

    def all(self):
        while not self._queue.empty():
            sol = self._search()
            self._solutions.append(sol)
            if len(self._solutions) >= self.option('max_solutions'):
                break
        return self._solutions

    def _first_round(self, start_n, start_ei):
        mlib_indexes = self._selector.get_connected_libs(-1)
        for ntype in mlib_indexes:
            motif_states = self._selector.get_motif_states(ntype)
            for ms in motif_states:
                new_ms = motif_state.Motif.copy(ms)
                motif_state.align_motif_state(start_n.state.get_end(start_ei), new_ms)
                self._test_node.update_ms(new_ms)
                score = self._scorer.score(self._test_node)
                new_node = MotifStateSearchNode(new_ms, start_n, start_ei, ntype, score)
                self._queue.put(new_node, score)

    def _search(self):
        accept_score, max_node_level, sterics, max_size, min_size = self._get_local_variables()
        best = 1000000000
        count = 0
        while not self._queue.empty():
            current = self._queue.get()
            score = self._scorer.accept_score(current)
            count += 1
            if score < best:
                best = score
                print best

            if count > 1000:
                break

            if score < accept_score:
                if current.size < min_size:
                    continue
                return MotifStateSearchSolution(current, score)

            if current.level+1 > max_node_level:
                continue

            mlib_indexes = self._selector.get_connected_libs(current.ntype)
            self._test_node.update_parent(current)
            for i, end in enumerate(current.state.iter_ends()):
                if i == 0:
                    continue
                for ntype in mlib_indexes:
                    motif_states = self._selector.get_motif_states(ntype)
                    for ms in motif_states:
                        motif_state.align_motif_state(end, ms)
                        self._test_node.update_ms(ms)
                        score = self._scorer.score(self._test_node)
                        #print score, current.score
                        if score > current.score:
                            continue
                        new_ms = motif_state.Motif.copy(ms)
                        new_node = MotifStateSearchNode(new_ms, current, i, ntype, score)
                        self._queue.put(new_node, score)

                        #print current.state.name, ms.name, score
        return None


class MotifStateSearchNode(object):
    def __init__(self, state, parent, parent_end_index, ntype, score):
        self.parent = parent
        self.parent_end_index = parent_end_index
        self.ss_score = state.score
        self.size =   state.num_res()
        self.ntype = ntype
        self.state = state
        if parent is None:
            self.level = 1
        else:
            self.level = self.parent.level + 1
            self.ss_score += self.parent.ss_score
            self.size += self.parent.size
        self.score = score

    def update_parent(self, parent):
        self.parent = parent
        self.ss_score = parent.ss_score
        self.size = parent.size

    def update_ms(self, ms):
        self.state = ms
        if self.parent is None:
            return
        self.level = self.parent.level + 1
        self.ss_score = self.state.score + self.parent.ss_score
        self.size = self.state.num_res() + self.parent.size


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

    def to_mst(self, rm):
        mst = motif_state_tree.MotifStateTree(rm)
        mst.set_sterics(0)
        for i, n in enumerate(self.path):
            new_ms = motif_state.Motif.copy(n.state, new_uuid=1)
            if i == 0:
                mst.add_state(new_ms)
            else:
                j = mst.add_state(new_ms, -1, n.parent_end_index)
                if j == -1:
                    raise ValueError("something went horribly wrong, cannot build solution")
        return mst

    def to_motif_tree(self, rm):
        return self.to_mst(rm).to_motif_tree()

