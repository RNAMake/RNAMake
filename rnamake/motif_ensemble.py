import random
import os
import settings
import motif

class MotifEnsembleMember(object):
    __slots__ = [
        '_motif',
        '_energy'
    ]

    def __init__(self, motif, energy):
        self._motif, self._energy = motif, energy

    def to_str(self):
        return self.motif.to_str() + "#" + str(self.energy)

    def copy(self):
        member_copy = MotifEnsembleMember(self.motif.copy(), self.energy)
        return member_copy

    @property
    def motif(self):
        return self._motif

    @property
    def energy(self):
        return self._energy


class MotifEnsemble(object):
    __slots__ = [
        "_end_id",
        "_members",
        "_block_end_add"
    ]

    def __init__(self, motifs, energies):
        if len(motifs) == 0:
            raise ValueError("must supply atleast one motif")

        if len(motifs) != len(energies):
            raise ValueError("must supply the same number of motifs and energies")

        self._members = []
        for i, m in enumerate(motifs):
            ms = MotifEnsembleMember(m, energies[i])
            self._members.append(ms)

        self._members.sort(key=lambda x: x.energy, reverse=False)
        self._block_end_add = self._members[0].motif.block_end_add
        self._end_id = self._members[0].motif.get_end_id(0)

    def __len__(self):
        return len(self._members)

    def __iter__(self):
        return self._members.__iter__()

    def copy(self):
        me_copy = MotifEnsemble()
        me_copy.id = self.id
        me_copy.block_end_add = self.block_end_add
        members = [mem.copy() for mem in self.members]
        me_copy.members = self.members

        return me_copy

    def most_populated_state(self):
        pass

    def get_random_member(self):
        return random.choice(self.members)

    def to_str(self):
        s = self._end_id + "{" + str(self._block_end_add) + "{"
        for ms in self._members:
            s += ms.to_str() + "{"
        return s

    def to_file(self, name="test.me"):
        s = self.to_str()
        f = open(name, "w")
        f.write(s)
        f.close()

    def get_state(self):
        mse = MotifStateEnsemble()
        motif_states = []
        energies = []

        for mem in self.members:
            motif_states.append(mem.motif.get_state())
            energies.append(mem.energy)

        mse.setup(self.id, motif_states, energies)
        return mse

    def to_pdb(self, name="test.pdb"):
        f = open (name, "w")
        for i, mem in enumerate(self.members):
            f.write("MODEL " + str(i+1) + "\n")
            f.write(mem.motif.to_pdb_str())
            f.write("ENDMDL\n")
        f.close()

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
        self.block_end_add = 0

    def setup(self, id, motif_states, energies):
        self.id = id
        for i in range(len(motif_states)):
            self.members.append((MotifStateEnsembleMember(motif_states[i], energies[i])))
        self.members.sort(key = lambda x : x.energy, reverse=False)
        self.block_end_add = self.members[0].motif_state.block_end_add

    def copy(self):
        mes_copy = MotifStateEnsemble()

        mes_copy.id =self.id
        mes_copy.block_end_add = self.block_end_add
        members = []
        for mem in self.members:
            members.append(mem.copy())
        mes_copy.members = members

        return mes_copy

    def to_str(self):
        s = self.id + "{" + str(self.block_end_add) + "{"
        for ms in self.members:
            s += ms.to_str() + "{"
        return s

    def get_random_member(self):
        return random.choice(self.members)

    def update_res_uuids(self, res):
        for mem in self.members:
            mem.motif_state.update_res_uuids(res)


def file_to_motif_ensemble(path):
    f = open(path)
    s = f.readline()
    f.close()

    return str_to_motif_ensemble(s)


def str_to_motif_ensemble(s):
    me = MotifEnsemble()
    spl = s.split("{")
    members = []
    me.id = spl.pop(0)
    me.block_end_add = int(spl.pop(0))
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
    spl = s.split("{")
    members = []
    mes.id = spl.pop(0)
    mes.block_end_add = int(spl.pop(0))
    for s in spl[:-1]:
        spl2 = s.split("#")
        m = motif.str_to_motif_state(spl2[0])
        energy = float(spl2[1])
        ms = MotifStateEnsembleMember(m, energy)
        members.append(ms)
    mes.members = members
    return mes


def motif_state_to_motif_state_ensemble(ms):
    mse = MotifStateEnsemble()
    mse.setup(ms.end_ids[0], [ms], [1])
    return mse

