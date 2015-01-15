import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Blue Paired with Yellow"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_blue_paired_with_yellow"
		self.default_params_ = []
		self.code_length_ = 10
		self.publishable_ = True


	def score(self, design, params):
		sequence = design['sequence']
		secstruct = design['secstruct']

		score = 100

		for ii in range(2,len(secstruct)- 20):
			if(secstruct[ii] == "."):
				if(sequence[ii] == "U"):
					score -= 1

		return score

