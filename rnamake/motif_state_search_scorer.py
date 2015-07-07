import util
import math


class MotifTreeStateSearchScorer(object):
    def __init__(self, target=None):
        self.target = None
        self.target_flip = None
        if target is not None:
            self.set_target(target)


    def set_target(self, target):
        self.target = target
        self.target_flip = target.copy()
        self.target_flip.flip()

    def score(self, node):
        raise ValueError("cannot call base MotifTreeStateSearchScorer")

    def accept_score(self, node):
        best_score = 10000
        for i, bp_state in enumerate(node.cur_state.end_states):
            if i == 0:
                continue
            score = util.distance(bp_state.d, self.target.d)

            r_diff       = util.matrix_distance(bp_state.r, self.target.r)
            r_diff_flip  = util.matrix_distance(bp_state.r, self.target_flip.r)

            if r_diff > r_diff_flip:
                r_diff = r_diff_flip

            score += 2*r_diff
            if score > best_score:
                best_score = score
        return score


class MTSS_GreedyBestFirstSearch(MotifTreeStateSearchScorer):
    def __init__(self, target=None):
        super(self.__class__, self).__init__(target)

    def score(self, node):
        best_score = 1000

        for i, bp_state in enumerate(node.cur_state.end_states):
            if i == 0:
                continue

            score = new_score_function(bp_state, self.target,
                                        self.target_flip)

            if score < best_score:
                best_score = score

        return best_score





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

