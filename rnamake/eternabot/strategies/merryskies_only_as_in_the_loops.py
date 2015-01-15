import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Only As in the loops"
		self.author_ = "merryskies"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_beta_test_experiment_title_only_as_in_the_loops"
		self.default_params_ = []
		self.code_length_ = 7
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False

	def score(self, design, params):
		score = 100
		secstruct = design['secstruct']
		sequence = design['sequence']
		for ii in range(0,len(secstruct)):
			if(secstruct[ii] == "." and sequence[ii] != "A"):
				score -= 1
		return score
