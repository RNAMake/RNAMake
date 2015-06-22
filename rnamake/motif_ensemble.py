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

    def get_most_populated_state(self):
        pass

    def get_random_state(self):
        pass

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
