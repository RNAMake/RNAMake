from . import option
from . import base
from . import motif_tree

class MotifTreePrecomputer(base.Base):
    def __init__(self):
        self.is_setup, self.mlib = 0, None
        self.setup_options_and_constraints()

    def setup_options_and_constraints(self):
        options = { 'data_output'     : 'text',
                    'max_bps_per_end' : 11,
                    'min_bps_per_end' : 0,
                    'clash_radius'    : motif_tree.MotifTree().clash_radius
                    'flip'            : None,
                    'motif_pos'       : None }

        self.options = option.Options(options)
        self.constraints = {}

    def precompute_library(self, mlib, name=None, **options):
        self.options.dict_set(options)
        for m in mlib.motifs():
            pass

    def precompute_motif(self):
        pass




