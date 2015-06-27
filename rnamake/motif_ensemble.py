import random
import os
import settings
import motif

class MotifEnsembleMember(object):
    __slots__ = ['motif', 'energy']

    def __init__(self, motif, energy):
        self.motif, self.energy = motif, energy

    def to_str(self):
        return self.motif.to_str() + "#" + str(self.energy)

    def copy(self):
        member_copy = MotifEnsembleMember(self.motif.copy(), self.energy)
        return member_copy


class MotifEnsemble(object):
    def __init__(self):
        self.id = ""
        self.members = []

    def setup(self, id, motifs, energies):
        self.id = id
        self.members = []
        for i, m in enumerate(motifs):
            ms = MotifEnsembleMember(m, energies[i])
            self.members.append(ms)

        self.members.sort(key = lambda x : x.energy, reverse=False)

    def copy(self):
        me_copy = MotifEnsemble()
        me_copy.id = self.id
        members = [mem.copy() for mem in self.members]
        me_copy.members = self.members

        return me_copy

    def most_populated_state(self):
        pass

    def get_random_state(self):
        pass

    def to_str(self):
        s = self.id + "$"
        for ms in self.members:
            s += ms.to_str() + "$"
        return s

    def get_state(self):
        mse = MotifStateEnsemble()
        motif_states = []
        energies = []

        for mem in self.members:
            motif_states.append(mem.motif.get_state())
            energies.append(mem.energy)

        mse.setup(self.id, motif_states, energies)
        return mse


class MotifStateEnsembleMember(object):
    __slots__ = ['motif_state', 'energy']

    def __init__(self, motif_state, energy):
        self.motif_state, self.energy = motif_state, energy

    def to_str(self):
        return self.motif_state.to_str() + "#" + str(self.energy)

    def copy(self):
        member_copy = MotifStateEnsembleMember(self.motif_state.copy(), self.energy)
        return member_copy


class MotifStateEnsemble(object):
    def __init__(self):
        self.id = ""
        self.members = []

    def setup(self, id, motif_states, energies):
        self.id = id
        for i in range(len(motif_states)):
            self.members.append((MotifStateEnsembleMember(motif_states[i], energies[i])))
        self.members.sort(key = lambda x : x.energy, reverse=False)

    def copy(self):
        mes_copy = MotifStateEnsemble()

        mes_copy.id =self.id
        members = []
        for mem in self.members:
            members.append(mem.copy())
        mes_copy.members = members

        return mes_copy

    def to_str(self):
        s = self.id + "$"
        for ms in self.members:
            s += ms.to_str() + "$"
        return s


def str_to_motif_ensemble(s):
    me = MotifEnsemble()
    spl = s.split("$")
    members = []
    me.id = spl.pop(0)
    for s in spl[:-1]:
        spl2 = s.split("#")
        m = motif.str_to_motif(spl2[0])
        energy = float(spl2[1])
        ms = MotifEnsembleMember(m, energy)
        members.append(ms)
    me.members = members
    return me

def str_to_motif_state_ensemble(s):
    mes = MotifStateEnsemble()
    spl = s.split("$")
    members = []
    mes.id = spl.pop(0)
    for s in spl[:-1]:
        spl2 = s.split("#")
        m = motif.str_to_motif_state(spl2[0])
        energy = float(spl2[1])
        ms = MotifStateEnsembleMember(m, energy)
        members.append(ms)
    mes.members = members
    return mes
