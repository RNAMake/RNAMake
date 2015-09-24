import motif_state_search
import steric_lookup
import resource_manager as rm
import base
import option

class MotifPathBuilder (base.Base):
    def __init__(self):
        self.search = motif_state_search.MotifStateSearch()
        self.setup_options_and_constraints()

    def setup_options_and_constraints(self):
        search_options = self.search.options.valid_options()
        self.option_location = {}

        options = {}
        for o in search_options:
            v = self.search.option(o)
            options[o] = v
            self.option_location[o] = 'search'

        self.options = option.Options(options)

    def _set_search_options(self):
        for k,v in self.options.options.iteritems():
            if k in self.option_location:
                if self.option_location[k] == 'search':
                    self.search.option(k,v)

    def build(self, start, end, p=None):

        self._set_search_options()
        if p is not None:
            beads = p.beads
            sl = steric_lookup.StericLookup()
            sl.add_beads(beads)
            self.search.lookup = sl

        while not self.search.finished():
            solution = self.search(one=1)






