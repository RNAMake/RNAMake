from rnamake.eternabot import sequence_designer
import resource_manager as rm
import motif_state_tree, monte_carlo, motif_state_graph, motif_state_search_scorer, vienna
import basepair, util
from eternabot import sequence_designer
import base
import option

from collections import defaultdict
import random


class OptimizedSequence(object):
    def __init__(self, sequence, dist_score, eterna_score):
        self.sequence = sequence
        self.dist_score = dist_score
        self.eterna_score = eterna_score


class SequenceOptimizerScorer(object):
    def __init__(self):
        pass

    def score(self, msg, ss):
        pass


class ExternalTargetScorer(SequenceOptimizerScorer):
    def __init__(self, target, ni, ei):
        self.target = target
        self.ni = ni
        self.ei = ei

    def score(self, msg, ss):
        state = msg.get_node(self.ni).data.cur_state.end_states[self.ei]
        return self.target.diff(state)


class InternalTargetScorer(SequenceOptimizerScorer):
    def __init__(self, ni1, ei1, ni2, ei2):
        self.ni1 = ni1
        self.ni2 = ni2
        self.ei1 = ei1
        self.ei2 = ei2

    def score(self, msg, ss):
        state1 = msg.get_node(self.ni1).data.cur_state.end_states[self.ei1]
        state2 = msg.get_node(self.ni2).data.cur_state.end_states[self.ei2]

        return state1.diff(state2)


class ViennaFoldScorer(SequenceOptimizerScorer):
    def __init__(self, target_ss, start=None, end=None):
        self.v = vienna.Vienna()
        self.target_ss = target_ss
        self.start = start
        self.end = end

    def score(self, msg, ss):
        if self.start is None:
            r = self.v.fold(ss.sequence())
        else:
            r = self.v.fold(ss.sequence()[self.start:self.end])
        diff = 0
        for i in range(len(self.target_ss)):
            if self.target_ss[i] != r.structure[i]:
                diff += 1
        return diff

class ViennaCoFoldScorer(SequenceOptimizerScorer):
    def __init__(self, target_ss, start_1=None, end_1=None, start_2=None, end_2=None):
        self.v = vienna.Vienna()
        self.target_ss = target_ss
        self.start_1 = start_1
        self.end_1 = end_1
        self.start_2 = start_2
        self.end_2 = end_2

    def score(self, msg, ss):
        if self.start_1 is None:
            seq1 = ss.sequence()
            seq2 = seq1
        else:
            seq1 = ss.sequence()[self.start_1:self.end_1]
            seq2 = ss.sequence()[self.start_2:self.end_2]

        r = self.v.cofold(seq1+"&"+seq2)
        diff = 0
        for i in range(len(self.target_ss)):
            if self.target_ss[i] != r.structure[i]:
                diff += 1

        # check homodimers
        r1 = self.v.cofold(seq1+"&"+seq1)
        r2 = self.v.cofold(seq2+"&"+seq2)

        if r1.energy < r.energy:
            diff += abs(r1.energy - r.energy)*5
        if r2.energy < r.energy:
            diff += abs(r2.energy - r.energy)*5

        return diff


class MultiScorer(SequenceOptimizerScorer):
    def __init__(self):
        self.scorers = []
        self.weights = []

    def add_scorer(self, scorer, weight):
        self.scorers.append(scorer)
        self.weights.append(weight)

    def score(self, msg, ss):
        total = 0
        for i, s in enumerate(self.scorers):
            total += s.score(msg, ss)*self.weights[i]
        return total


class SequenceOptimizer3D(base.Base):

    class _DesignableBP(object):
        def __init__(self, bp):
            self.bp = bp
            self.last_state = ""
            self.m_id_bot = None
            self.m_id_top = None

        def update_state(self, bp_name):
            self.last_state = self.bp.res1.name+self.bp.res2.name
            self.bp.res1.name = bp_name[0]
            self.bp.res2.name = bp_name[1]

        def revert_state(self):
            self.bp.res1.name = self.last_state[0]
            self.bp.res2.name = self.last_state[1]

        def name(self):
            return self.bp.res1.name+self.bp.res2.name

    def __init__(self, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.sequence_designer = sequence_designer.SequenceDesigner()
        self.scorer = None

        self.target_state = None
        self.target_state_flip = None
        self.possible_bps = ["AU", "UA", "GC", "CG"]

        self.steps = 0

    def setup_options_and_constraints(self):
        options = { 'cutoff'        : 5,
                    'solutions'     : 1,
                    'eterna_cutoff' : -1,
                    'return_lowest' : 0,
                    'verbose'       : 0,
                    'max_steps'     : 500}

        self.options = option.Options(options)

    def set_scorer(self, scorer):
        self.scorer = scorer

    def _update_designable_bp(self, d_bp, msg, ss):
        if d_bp.m_id_bot is not None:
            ss.update_motif(d_bp.m_id_bot)
            m = ss.motif(d_bp.m_id_bot)
            n = msg.get_node(uuid=d_bp.m_id_bot)
            msg.replace_state(n.index, rm.manager.get_bp_step_state(m.end_ids[0]))

        if d_bp.m_id_top is not None:
            ss.update_motif(d_bp.m_id_top)
            m = ss.motif(d_bp.m_id_top)
            n = msg.get_node(uuid=d_bp.m_id_top)
            msg.replace_state(n.index, rm.manager.get_bp_step_state(m.end_ids[0]))

    def _get_designable_bps(self, ss):
        designable_bps = []
        for bp in ss.basepairs:
            bp_name = bp.res1.name+bp.res2.name
            if bp_name == "NN":
                designable_bps.append(self._DesignableBP(bp))

        for d_bp in designable_bps:
            for m in ss.motifs:
                if m.name != "HELIX.IDEAL":
                    continue
                if m.ends[0] == d_bp.bp:
                    d_bp.m_id_bot = m.id
                if m.ends[1] == d_bp.bp:
                    d_bp.m_id_top = m.id
            state = random.choice(self.possible_bps)
            d_bp.bp.res1.name = state[0]
            d_bp.bp.res2.name = state[1]

        return designable_bps

    def _validate_sequence(self, msg, ss):
        s1 = msg.to_motif_graph().secondary_structure().sequence()
        s2 = ss.sequence()

        for j in range(len(s2)):
            if s1[j] != s2[j]:
                print s1
                print s2
                raise ValueError(
                    "sequences are out of sync: something went really wrong in sequence "
                    "optimization")
        return s1

    def _initiate_sequence_in_msg(self, msg, ss):
        for m in ss.motifs:
            ss.update_motif(m.id)

        for m in ss.motifs:
            n = msg.get_node(uuid=m.id)

            if n.data.name() != "HELIX.IDEAL":
                continue

            msg.replace_state(n.index, rm.manager.get_bp_step_state(m.end_ids[0]))

    def get_optimized_sequences(self, mg, scorer=None):
        if scorer is not None:
            self.set_scorer(scorer)

        if self.scorer is None:
            raise ValueError(
                "cannot run get_optimized_sequences without scorer, either supply here "
                "or use set_scorer")

        org_ss = mg.secondary_structure()
        org_score = 0
        ss = mg.designable_secondary_structure()

        if len(ss.chains()) > 1:
            raise ValueError(
                "cannot perform sequence optmization with more than one chain")

        designable_bps = self._get_designable_bps(ss)

        msg = motif_state_graph.MotifStateGraph(mg)
        self._initiate_sequence_in_msg(msg, ss)
        last_score = self.scorer.score(msg, ss)

        mc = monte_carlo.MonteCarlo()
        mc.temperature = 1

        best_states = []
        best_seq = ""
        best = 1000
        cutoff = self.option('cutoff')
        solutions = []
        for i in range(1, self.option('max_steps')):
            d_bp = random.choice(designable_bps)

            new_bp_state = random.choice(self.possible_bps)
            d_bp.update_state(new_bp_state)

            self._update_designable_bp(d_bp, msg, ss)

            new_score = self.scorer.score(msg, ss)

            if mc.accept(last_score, new_score):
                last_score = new_score
            else:
                d_bp.revert_state()
                self._update_designable_bp(d_bp, msg, ss)
                continue

            if best > new_score:
                if self.option('verbose'):
                    print best, new_score
                best_seq = ss.sequence()
                best = new_score

            if cutoff < new_score:
                continue

            eterna_score = self.sequence_designer.design(ss.dot_bracket(),
                                                         ss.sequence())[0].score

            if eterna_score > org_score:
                s = self._validate_sequence(msg, ss)
                solutions.append(
                    OptimizedSequence(s, new_score, eterna_score))

                if len(solutions) >= self.option('solutions'):
                    return solutions


        if len(solutions) == 0 and self.option('return_lowest'):
            return [OptimizedSequence(best_seq, best, -1)]


        return solutions

    def get_optimized_mg(self, mg, scorer=None):
        if scorer is not None:
            self.set_scorer(scorer)

        if self.scorer is None:
            raise ValueError(
                "cannot run get_optimized_sequences without scorer, either supply here "
                "or use set_scorer")

        org_ss = mg.secondary_structure()
        org_score = 0
        ss = mg.designable_secondary_structure()

        if len(ss.chains()) > 1:
            raise ValueError(
                "cannot perform sequence optmization with more than one chain")

        designable_bps = self._get_designable_bps(ss)

        msg = motif_state_graph.MotifStateGraph(mg)
        self._initiate_sequence_in_msg(msg, ss)
        last_score = self.scorer.score(msg)

        mc = monte_carlo.MonteCarlo()
        mc.temperature = 1

        best_states = []
        best_seq = ""
        best = 1000
        cutoff = self.option('cutoff')
        solutions = []
        for i in range(1, self.option('max_steps')):
            d_bp = random.choice(designable_bps)

            new_bp_state = random.choice(self.possible_bps)
            d_bp.update_state(new_bp_state)

            self._update_designable_bp(d_bp, msg, ss)

            new_score = self.scorer.score(msg)

            if mc.accept(last_score, new_score):
                last_score = new_score
            else:
                d_bp.revert_state()
                self._update_designable_bp(d_bp, msg, ss)
                continue

            if best > new_score:
                if self.option('verbose'):
                    print best, new_score
                best_seq = ss.sequence()
                best = new_score

            if cutoff < new_score:
                continue

            eterna_score = self.sequence_designer.design(ss.dot_bracket(),
                                                         ss.sequence())[0].score

            if eterna_score > org_score:
                s = self._validate_sequence(msg, ss)
                return msg.to_motif_graph()


        if len(solutions) == 0 and self.option('return_lowest'):
            return msg.to_motif_graph()