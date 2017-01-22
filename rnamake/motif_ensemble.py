import motif
import motif_state
from primitives.ensemble import Ensemble, EnsembleMember


class MotifEnsemble(Ensemble):
    __slots__ = [
        "_end_id",
        "_members",
        "_block_end_add"
    ]

    def __init__(self, motifs, energies):
        if len(motifs) != len(energies):
            raise ValueError("must supply the same number of motifs and energies")

        members = []
        for i, m in enumerate(motifs):
            ms = EnsembleMember(m, energies[i])
            members.append(ms)

        members.sort(key=lambda x: x.energy, reverse=False)
        super(self.__class__, self).__init__(members)

    @classmethod
    def from_str(cls, s, rts):
        spl = s.split("{")
        motifs = []
        energies = []
        for s in spl[:-1]:
            spl2 = s.split("#")
            m = motif.Motif.from_str(spl2[0], rts)
            energy = float(spl2[1])
            motifs.append(m)
            energies.append(energy)
        return cls(motifs, energies)

    @classmethod
    def copy(cls, me):
        motifs = []
        energies = []
        for mem in me:
            motifs.append(motif.Motif.copy(mem.motif))
            energies.append(mem.energy)
        return cls(motifs, energies)

    def to_str(self):
        s = ""
        for ms in self._members:
            s += ms.to_str() + "{"
        return s

    def to_file(self, name="test.me"):
        s = self.to_str()
        f = open(name, "w")
        f.write(s)
        f.close()

    def get_state(self):
        motif_states = []
        energies = []

        for mem in self._members:
            motif_states.append(mem.motif.get_state())
            energies.append(mem.energy)

        return motif_state.MotifEnsemble(motif_states, energies)

    def to_pdb(self, name="test.pdb"):
        f = open (name, "w")
        for i, mem in enumerate(self._members):
            f.write("MODEL " + str(i+1) + "\n")
            f.write(mem.motif.to_pdb_str())
            f.write("ENDMDL\n")
        f.close()


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



