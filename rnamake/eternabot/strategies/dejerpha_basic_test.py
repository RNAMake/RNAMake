import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
    def __init__(self):
        strategy_template.Strategy.__init__(self)

        self.title_ = "A basic test"
        self.author_ = "dejerpha"
        self.url_ = "http://getsatisfaction.com/eternagame/topics/a_basic_test"
        self.default_params_ = [0.5, 100, -1.5, 1, 1, 77, 97]
        self.code_length_ = 20
        self.publishable_ = True
        #self.denormalized_ = True
        self.comprehensive_ = True

    def score(self, design, params):
        sequence = design['sequence']
        seqlen = len(sequence)
        seq_range = range(0,len(sequence))
        fe = design['fe']
        mp = design['meltpoint']

        total_pairs = design['gu'] + design['gc'] + design['ua']

        score = 100
        if(total_pairs > 0):
            score -= abs(float(design['ua'])/total_pairs - params[0]) * params[1]

        target_fe = params[2] * total_pairs

        score -= abs(target_fe - fe) * params[3]

        if(mp < params[5]):
            score -= abs(mp -params[5]) * params[4]
        elif(mp > params[6]):
            score -= abs(mp - params[6]) * params[4]


        return score

