import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
    def __init__(self):
        strategy_template.Strategy.__init__(self)

        self.title_ = "Clean plot, stack caps, and safe GC"
        self.author_ = "xmbrst"
        self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_clean_plot_stack_caps_and_safe_gc"
        self.default_params_ = [1.0, 1.0, 2.0, 0.8]
        self.code_length_ = 70
        self.publishable_ = True
        self.denormalized_ = True
        self.comprehensive_ = True

    def score(self, design, params):

        penalty = 0.0

        n = len(design['sequence'])
        sequence = design['sequence']
        npairs = design['gc'] + design['gu'] + design['ua']
        dotplot = design['dotplot']
        pairmap = design['pairmap']
        elements = design['secstruct_elements']

        for ii in range(0,len(dotplot)):
            i_index = dotplot[ii][0]
            j_index = dotplot[ii][1]

            if(pairmap[i_index] != j_index):
                penalty += dotplot[ii][2]

        if(npairs > 0):
            plotscore = (1.0 - (penalty / npairs))
        else:
            plotscore = 0
        gc_penalty = 0
        if(npairs > 0):
            if(float(design['gc'])/npairs > params[3]):
                gc_penalty = 1

        cap_score = 0
        stack_count = 0
        for ii in range(0,len(elements)):
            if(elements[ii].type_ == eterna_utils.RNAELEMENT_STACK):
                stacklen = len(elements[ii].indices_)/2
                stack_count += 1
                if(stacklen == 1):
                    pair = elements[ii].get_pair_from_stack(0, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 1
                elif(stacklen == 2):
                    pair = elements[ii].get_pair_from_stack(0, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 0.5

                    pair = elements[ii].get_pair_from_stack(1, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 0.5
                elif(stacklen == 3):

                    pair = elements[ii].get_pair_from_stack(0, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 0.4

                    pair = elements[ii].get_pair_from_stack(1, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 0.4

                    pair = elements[ii].get_pair_from_stack(2, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 0.4

                else:
                    pair = elements[ii].get_pair_from_stack(0, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 1 / 3.0

                    pair = elements[ii].get_pair_from_stack(1, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 1 / 6.0

                    pair = elements[ii].get_pair_from_stack(stacklen-2, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 1 / 6.0

                    pair = elements[ii].get_pair_from_stack(stacklen-1, sequence)
                    if(pair == "GC" or pair == "CG"):
                        cap_score += 1 / 3.0


        if(stack_count > 0):
            cap_score = float(cap_score)/stack_count

        return (2.0 + cap_score * params[1] + plotscore * params[0] - gc_penalty * params[2]) * 25;
