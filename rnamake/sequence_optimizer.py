from rnamake.eternabot import sequence_designer
import resource_manager as rm
import motif_state_tree, monte_carlo
from rnamake import basepair
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


class SequenceOptimizer3D(base.Base):
    def __init__(self, **options):
        self.opt_nodes = []
        self.bps = {}
        self.seqs = []
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.sequence_designer = sequence_designer.SequenceDesigner()

    def setup_options_and_constraints(self):
        options = { 'cutoff'       : 5,
                    'solutions'    : 10,
                    'eterna_cutoff' : -1}

        self.options = option.Options(options)

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


    def _update_designable_bp(self, d_bp, mst, new_bp_state):
        d_bp.state = new_bp_state

        if d_bp.m_id_bot is not None:
            n = mst.get_node(uuid=d_bp.m_id_bot)
            spl = n.data.cur_state.name.split("=")
            new_name = new_bp_state + "=" + spl[1]
            mst.replace_state(n.index, rm.manager.get_state(name=new_name))

        if d_bp.m_id_top is not None:
            n = mst.get_node(uuid=d_bp.m_id_top)
            spl = n.data.cur_state.name.split("=")
            new_name = spl[0] + "=" + new_bp_state
            mst.replace_state(n.index, rm.manager.get_state(name=new_name))


    def get_optimized_sequences(self, mt, target_bp, ni, ei):
        org_ss = mt.secondary_structure()
        org_score = self.sequence_designer.design(org_ss.dot_bracket(), org_ss.sequence())[0].score

        ss = mt.designable_secondary_structure()
        designable_bps = []
        for bp in ss.basepairs:
            bp_name = bp.res1.name+bp.res2.name
            if bp_name == "NN":
                designable_bps.append(self._DesignableBP(bp))

        possible_bps = ["AU", "UA", "GC", "CG"]

        for d_bp in designable_bps:
            for m in ss.motifs:
                if m.name != "HELIX.IDEAL":
                    continue
                if m.ends[0] == d_bp.bp:
                    d_bp.m_id_bot = m.id
                if m.ends[1] == d_bp.bp:
                    d_bp.m_id_top = m.id
            state = random.choice(possible_bps)
            d_bp.bp.res1.name = state[0]
            d_bp.bp.res2.name = state[1]

        mst = motif_state_tree.MotifStateTree(mt=mt)
        for m in ss.motifs:
            if m.name != "HELIX.IDEAL":
                continue

            n = mst.get_node(uuid=m.id)
            name =  m.ends[0].res1.name+m.ends[0].res2.name+"="
            name += m.ends[1].res1.name+m.ends[1].res2.name
            mst.replace_state(n.index, rm.manager.get_state(name=name))

        target_state = target_bp.state()
        target_state_flip = target_bp.state()
        target_state_flip.flip()

        last_score = target_state.diff(mst.get_node(ni).data.cur_state.end_states[ei])

        mc = monte_carlo.MonteCarlo()
        mc.temperature = 5

        best_states = []
        best = 1000
        cutoff = self.option('cutoff')
        solutions = []
        for i in range(1,5000):
            if i % 500 == 0:
                print i, best
                for j, s in enumerate(best_states):
                    name = s[1]
                    if name[2] == "=":
                        mst.replace_state(s[0], rm.manager.get_state(name=name))
                mc.temperature *= 0.8

            d_bp = random.choice(designable_bps)

            new_bp_state = random.choice(possible_bps)
            d_bp.update_state(new_bp_state)

            self._update_designable_bp(d_bp, mst, new_bp_state)

            new_score = target_state.diff(mst.get_node(ni).data.cur_state.end_states[ei])

            #print last_score, new_score
            if mc.accept(last_score, new_score):
                last_score = new_score
            else:
                self._update_designable_bp(d_bp, mst, d_bp.last_state)
                continue


            if best > new_score:
                print best
                best = new_score
                best_states = []
                for n in mst:
                    best_states.append([n.index,
                                        n.data.cur_state.name,
                                        n.data.cur_state.end_names[0]])


            if cutoff < new_score:
                continue

            eterna_score = self.sequence_designer.design(ss.dot_bracket(),
                                                         ss.sequence())[0].score
            if eterna_score > org_score:
                solutions.append(
                    OptimizedSequence(ss.sequence(), new_score, eterna_score))

                if len(solutions) > self.option('solutions'):
                    return solutions



        return solutions




