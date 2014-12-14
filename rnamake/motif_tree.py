from . import base
from . import option
from . import settings
from . import motif

class MotifTree(base.Base):
    def __init__(self, motif=None, **options):
        self.setup_options_and_constraints()
        self.options.dict_add(options)

        if motif is None:
            head = None

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1,
                    'full_beads_first_res' : 1,
                    'ideal_bp_score'       : -2 # TODO fix
                    }

        self.options = option.Options(options)
        self.constraints = {}

    def _get_default_head(self):
        mdir = settings.RESOURCES_PATH + "/start"
        m = motif.Motif(mdir)
        return MotifTreeNode(m, 0)

def MotifTreeNode(object):
    def __int__(self, motif, level):
        pass
