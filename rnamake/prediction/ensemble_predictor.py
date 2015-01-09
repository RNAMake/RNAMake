import motif_ensemble_tree
import rnamake.secondary_structure_tree as secondary_structure_tree
import rnamake.motif_tree_state as motif_tree_state

class EnsemblePredictor(object):
    def __init__(self, sequence=None, structure=None, met=None):
        self.met = motif_ensemble_tree.MotifEnsembleTree()
        if met is not None:
            self.met = met
        self.mtst = motif_tree_state.MotifTreeStateTree()


