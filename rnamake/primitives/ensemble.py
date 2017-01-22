import random

class Ensemble(object):
    __slots__ = [
        '_members',
        '_end_id',
        '_block_end_add'
    ]

    def __init__(self, members):
        if len(members) == 0:
            raise ValueError("cannot have an ensemble with 0 members")

        self._members = members
        self._end_id  = self._members[0].motif.get_end_id(0)
        self._block_end_add = self._members[0].motif.block_end_add

    def __len__(self):
        return len(self._members)

    def __iter__(self):
        return self._members.__iter__()

    def get_random_member(self):
        return random.choice(self._members)

    def get_member(self, i):
        return self._members[i]

    def get_energies(self):
        energies = []
        for mem in self._members:
            energies.append(mem.energy)
        return energies

    def get_motifs(self):
        motifs = []
        for mem in self._members:
            motifs.append(mem.motifs)
        return motifs

    @property
    def end_id(self):
        return self._end_id

    @property
    def block_end_add(self):
        return self._block_end_add


class EnsembleMember(object):
    __slots__ = [
        '_motif',
        '_energy'
    ]

    def __init__(self, motif, energy):
        self._motif, self._energy = motif, energy

    def to_str(self):
        return self._motif.to_str() + "#" + str(self._energy)

    @property
    def motif(self):
        return self._motif

    @property
    def energy(self):
        return self._energy
