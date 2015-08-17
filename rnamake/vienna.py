import settings
import os

class ViennaResults(object):
    def __init__(self, structure, energy, ensemble_prob, ensemble_diversity):
        self.structure, self.energy = structure, energy
        self.ensemble_prob, self.ensemble_diversity = ensemble_prob, ensemble_diversity

class Vienna(object):
    def __init__(self):
        self.bin_path = settings.VIENNA_BIN

    def fold(self, seq):
        if len(seq) == 0:
            raise ValueError("must supply a sequence longer then 0")
        os.system("echo \""+seq+"\" | "+self.bin_path+"RNAfold -p > rnafold_dump")
        f = open("rnafold_dump")
        lines = f.readlines()
        f.close()
        try:
            last_line = lines.pop()
        except:
            return None
        spl = last_line.split()
        ensemble_prob = float(spl[6][:-1])
        ensemble_diversity = float(spl[-1])
        spl = lines[1].split()
        spl2 = lines[1].split("(")
        structure = spl[0]
        energy = float(spl2[-1][:-2].rstrip())
        results = ViennaResults(structure, energy, ensemble_prob,
                                ensemble_diversity)
        os.remove("rnafold_dump")
        return results


    def cofold(self, seq):
        if len(seq) == 0:
            raise ValueError("must supply a sequence longer then 0")
        os.system("echo \""+seq+"\" | "+self.bin_path+"RNAcofold -p > rnafold_dump")
        f = open("rnafold_dump")
        lines = f.readlines()
        f.close()
        try:
            last_line = lines.pop()
        except:
            return None
        spl = last_line.split()
        ensemble_prob = float(spl[6][:-1])
        ensemble_diversity = float(spl[-1])
        spl = lines[1].split()
        spl2 = lines[1].split("(")
        structure = spl[0]
        energy = float(spl2[-1][:-2].rstrip())
        results = ViennaResults(structure, energy, ensemble_prob,
                                ensemble_diversity)
        os.remove("rnafold_dump")
        return results
