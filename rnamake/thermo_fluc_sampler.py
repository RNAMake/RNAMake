import base
import copy
import option
import random
import math
import basepair
import util
import motif_outputer

class ThermoFlucSampler(base.Base):
    def __init__(self, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.mset = None
        self.mst = None

    def setup_options_and_constraints(self):
        options = { 'temperature' : 298.15 }

        self.options = option.Options(options)
        self.constraints = {}

    def setup(self, mset):
        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        # kB.T at room temperature (25 degree Celsius)
        self.kBT = kB * self.option('temperature')

        self.mset = mset
        self.mst  = mset.to_mst()
        self.states = [ 0 for x in range(len(self.mst)) ]

    def sample(self, i):
        for j in range(i):
            self.next()

    def next(self):
        node_num  = random.randint(2, len(self.mset)-1)
        mset_node = self.mset.get_node(node_num)
        mst_node  = self.mst.get_node(node_num)
        pos       = self.states[node_num]

        energy    = mset_node.data.members[pos].energy
        new_mem   = mset_node.data.get_random_member()

        if new_mem.energy < energy:
            return self.update(node_num, new_mem)

        score = math.exp((energy - new_mem.energy) / self.kBT)
        dice_roll = random.random()
        if dice_roll < score:
            return self.update(node_num, new_mem)

        return 0

    def update(self, node_num, new_mem):
        mset_node = self.mset.get_node(node_num)
        self.last_state_pos = self.states[node_num]
        self.states[node_num] = mset_node.data.members.index(new_mem)
        self.last_state = self.mst.get_node(node_num).data.ref_state
        self.last_num = node_num
        self.mst.replace_state(node_num, new_mem.motif_state)
        return 1

    def undo(self):
        node_num = self.last_num
        state = self.last_state
        self.states[node_num] = self.last_state_pos
        self.mst.replace_state(node_num, state)

    def to_pdb(self, name="test.pdb"):
        self.mst.to_pdb(name)


class ThermoFlucRelax(base.Base):
    def __init__(self):
        self.sampler = ThermoFlucSampler(temperature=1000.0)
        self.max_steps = 10000
        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        self.kBT = kB * 100

    def run(self, mset, ni_1, ni_2, ei_1, ei_2):
        self.sampler.setup(mset)
        steps = 0

        end_state_1 = self.sampler.mst.get_node(ni_1).data.cur_state.end_states[ei_1]
        end_state_2 = self.sampler.mst.get_node(ni_2).data.cur_state.end_states[ei_2]
        cur_diff = end_state_1.diff(end_state_2)
        new_diff = 0
        self.best = self.sampler.mst.copy()
        best_diff = 0

        #print len(self.sampler.mst)
        #self.sampler.mst.add_connection(0, 16, "A149-A154")
        #outputer = motif_outputer.MotifOutputer()

        while steps < self.max_steps:

            #print steps
            r = self.sampler.next()

            """if steps % 10 == 0:
                try:
                    mst2 = self.sampler.mst.copy()
                    p = mst2.to_pose()
                    outputer.add_motif(p, 0)

                except:
                    pass"""
            steps += 1

            if steps % 1000 == 0:
                self.kBT *= 0.5

            if r == 0:
                continue

            end_state_1 = self.sampler.mst.get_node(ni_1).data.cur_state.end_states[ei_1]
            end_state_2 = self.sampler.mst.get_node(ni_2).data.cur_state.end_states[ei_2]
            new_diff = end_state_1.diff(end_state_2)

            #if steps % 10 == 0:
            #    print cur_diff, new_diff, best_diff, steps

            if new_diff > best_diff:
                best_diff = new_diff
                self.best = self.sampler.mst.copy()

            if new_diff > cur_diff:
                cur_diff = new_diff
                continue

            score = math.exp((cur_diff - new_diff) / self.kBT)
            dice_roll = random.random()
            if dice_roll < score:
                cur_diff = new_diff
                continue

            self.sampler.undo()
        #print end_state_1.diff(end_state_2)
        #outputer.to_pdb()


    def to_pdb(self, name="test.pdb"):
        self.best.to_pdb(name)

    def write_pdbs(self):
        self.best.write_pdbs()


class ThermoFlucScorer(object):
    def __init__(self):
        pass

class FrameScorer(ThermoFlucScorer):
    def __init__(self):
        pass

    def score(self, state_1, state_2):
        frame_score = util.distance(state_1.d, state_2.d)
        r_diff = util.matrix_distance(state_1.r, state_2.r)
        state_2.flip()
        r_diff_flip = util.matrix_distance(state_1.r, state_2.r)
        state_2.flip()
        if r_diff > r_diff_flip:
            r_diff = r_diff_flip
        frame_score += r_diff
        return frame_score


class RMSDScorer(ThermoFlucScorer):
    pass


class ThermoFlucSimulation(base.Base):
    def __init__(self, **options):
        self.setup_options_and_constraints()
        self.scorer = FrameScorer()
        self.sampler = ThermoFlucSampler()

    def setup_options_and_constraints(self):
        options   = { 'temperature' : 298.15,
                      'steps'       : 100000,
                      'record'      : 1}

        opt_types =  { 'temperature' : 'sampler'  }

        self.options = option.Options(options)
        self.opt_types = opt_types
        self.constraints = {}

    def setup(self, mset, ni1, ni2, ei1, ei2, **options):
        self.options.dict_set(options)

        self.sampler.setup(mset)
        self.n1 = self.sampler.mst.get_node(ni1)
        self.n2 = self.sampler.mst.get_node(ni2)
        self.ei1, self.ei2 = ei1, ei2

        end_state_1 = self.n1.data.cur_state.end_states[self.ei1]
        end_state_2 = self.n2.data.cur_state.end_states[self.ei2]

    def run(self):

        step = 0
        max_steps = self.option('steps')
        record = self.option('record')

        if record:
            f = open("results.txt", "w")

        count = 0
        while step < max_steps:
            self.sampler.next()

            end_state_1 = self.n1.data.cur_state.end_states[self.ei1]
            end_state_2 = self.n2.data.cur_state.end_states[self.ei2]

            score = self.scorer.score(end_state_1, end_state_2)

            if score < 5:
                count += 1
                #d_diff = util.distance(end_state_1.d, end_state_2.d)
                #r_diff = end_state_1._rot_diff(end_state_2)
                #f.write(str(d_diff) + " " + str(r_diff) + "\n")

            if record:
                d_diff = util.distance(end_state_1.d, end_state_2.d)
                r_diff = end_state_1._rot_diff(end_state_2)
                f.write(str(d_diff) + " " + str(r_diff) + "\n")

            step += 1

        if record:
            f.close()

        print count








