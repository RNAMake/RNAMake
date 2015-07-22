import base
import option
import random
import math

class ThermoFlucSampler(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()
        self.mset = None
        self.mst = None

        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        self.kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1}

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

    def update(self, node_num, new_mem):
        mset_node = self.mset.get_node(node_num)
        self.last_state_pos = self.states[node_num]
        self.states[node_num] = mset_node.data.members.index(new_mem)
        self.last_state = self.mst.get_node(node_num).data.ref_state
        self.last_num = node_num
        self.mst.replace_state(node_num, new_mem.motif_state)

    def undo(self):
        node_num = self.last_num
        state = self.last_state
        self.states[node_num] = self.last_state_pos
        self.mst.replace_state(node_num, state)

    def to_pdb(self, name="test.pdb"):
        self.mst.to_pdb(name)


class ThermoFlucRelax(base.Base):
    def __init__(self):
        self.sampler = ThermoFlucSampler()

    def setup(self, mset):
        pass