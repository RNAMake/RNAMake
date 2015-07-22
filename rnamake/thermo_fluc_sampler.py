import base
import option
import random
import math
import basepair

class ThermoFlucSampler(base.Base):
    def __init__(self, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.mset = None
        self.mst = None

        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        # kB.T at room temperature (25 degree Celsius)
        self.kBT = kB * self.option('temperature')

    def setup_options_and_constraints(self):
        options = { 'temperature' : 298.15 }

        self.options = option.Options(options)
        self.constraints = {}

    def setup(self, mset):
        self.mset = mset
        self.mst  = mset.to_mst()
        self.states = [ 0 for x in range(len(self.mst)) ]

    def sample(self, i):
        for j in range(i):
            self.next()

    def next(self):
        node_num  = random.randint(1, len(self.mset)-1)
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
        self.max_steps = 1000

    def run(self, mset, ni_1, ni_2, ei_1, ei_2):
        self.sampler.setup(mset)
        steps = 0

        end_state_1 = self.sampler.mst.get_node(ni_1).data.cur_state.end_states[ei_1]
        end_state_2 = self.sampler.mst.get_node(ni_2).data.cur_state.end_states[ei_2]
        cur_diff = end_state_1.diff(end_state_2)
        new_diff = 1000

        while steps < self.max_steps:
            r = self.sampler.next()
            if r == 0:
                continue

            end_state_1 = self.sampler.mst.get_node(ni_1).data.cur_state.end_states[ei_1]
            end_state_2 = self.sampler.mst.get_node(ni_2).data.cur_state.end_states[ei_2]
            new_diff = end_state_1.diff(end_state_2)
            if new_diff > cur_diff:
                self.sampler.undo()
            else:
                cur_diff = new_diff
            steps += 1
        print cur_diff


    def to_pdb(self, name="test.pdb"):
        self.sampler.to_pdb(name)

    def write_pdbs(self):
        self.sampler.mst.write_pdbs()
