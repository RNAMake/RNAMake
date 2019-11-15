from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
    def __init__(self):
        strategy_template.Strategy.__init__(self)

        self.title_ = "Direction of GC-pairs in multiloops + neckarea"
        self.author_ = "Eli Fisker"
        self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_direction_of_gc_pairs_in_multiloops_neckarea"
        self.default_params_ = [-1, -2] # Penalty for wrong GC, penalty for non-GC
        self.code_length_ = 30
        self.publishable_ = True
        self.comprehensive_ = False
        self.denormalized_ = True

    def score(self, design, params):
        elements = design['secstruct_elements']

        penalty = 0
        multiloop_found = False

        for e_index, element in enumerate(elements):
            if (element.type_ != RNAELEMENT_LOOP): continue
            closing_pairs = element.get_loop_closing_pairs(design['sequence'],design['pairmap'])
            if len(closing_pairs) < 3: continue
            multiloop_found = True
            for index, pair in enumerate(closing_pairs):
                #if (e_index == 1) and (index == 0):
                #    print("made it")
                #if (e_index == 1) and (index == 0): continue # loop closing pair connected with neck is allowed to turn any direction
                if index == 0:
                    if pair == "CG": continue
                    elif pair == "GC": penalty += params[0]
                    else:
                        penalty += params[1]
                else:
                    if pair == "GC": continue
                    elif pair == "CG": penalty += params[0]
                    else:
                        penalty += params[1]
        if not multiloop_found: return 100
        return 100 + penalty
