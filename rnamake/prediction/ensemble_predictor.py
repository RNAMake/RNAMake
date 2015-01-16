import motif_ensemble_tree
import motif_ensemble
import rnamake.secondary_structure_tree as secondary_structure_tree
import rnamake.motif_tree_state as motif_tree_state
import random
import math

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1# use start += 1 to find overlapping matches


class EnsemblePredictor(object):
    def __init__(self, sequence, structure, met=None, motifs=None):
        if met is not None:
            self.met = met
            return
        self.ss_tree = secondary_structure_tree.SecondaryStructureTree(structure,
                                                                       sequence)

        """seq = motifs[0].sequence()
        seq_chains = seq.split("+")
        saved = []
        pos = 0
        print seq_chains
        print list(find_all(sequence, seq_chains[0]))
        print sequence.find(seq_chains[0])
        print sequence.find(seq_chains[1])
        print sequence.find(seq_chains[2])
        return"""
        self.met = motif_ensemble_tree.MotifEnsembleTree()
        helical_steps = []
        for n in self.ss_tree.nodes:
            if n.ss_type == "Basepair":
                helical_steps.append(n.bp_type)
            if n.ss_type == "Bulge":
                self.match_bulge_to_motif(n, motifs)
                exit()

        step_names = self.get_step_motifs(helical_steps)
        for name in step_names:
            self.met.add_ensemble(motif_ensemble.MotifEnsemble(name, 0, 0))


    def match_bulge_to_motif(self, node, motifs):
        print node.x_seq, node.y_seq
        print motifs
        exit()


    def get_step_motifs(self, steps):
        motif_names = []
        for i in range(1,len(steps)):
            full_step = steps[i-1] + "=" + steps[i]
            motif_names.append(full_step)
        return motif_names

    def sample(self, steps=100, output_pdbs=0):
        mtst = self.met.get_mtst()
        temp = 0.04
        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)
        for i in range(steps):
            node_num = random.randint(1, len(self.met.nodes)-1)
            met_node = self.met.nodes[node_num]
            mts_node = mtst.nodes[node_num]
            cpop = met_node.motif_ensemble.get_state(mts_node.mts.name).population
            ms = met_node.get_random_state()
            if ms.population < cpop:
                mtst.replace_state(mts_node, ms.mts, 0)
                accept_count += 1
                if output_pdbs:
                    mtst.to_pdb("step."+str(i)+".pdb")
                continue
            score = math.exp((cpop - ms.population)/kBT)
            dice_roll = random.random()
            if dice_roll < score:
                accept_count += 1
                mtst.replace_state(mts_node, ms.mts, 0)
                if output_pdbs:
                    mtst.to_pdb("step."+str(i)+".pdb")





