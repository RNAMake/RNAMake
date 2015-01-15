import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):

		strategy_template.Strategy.__init__(self)

		self.title_ = "Twisted Basepairs"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_twisted_basepairs"
		self.default_params_ = []
		self.code_length_ = 15
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False

	def score(self, design, params):
		score = 100

		elements = design['secstruct_elements']
		pairmap = design['pairmap']
		sequence = design['sequence']

		for ii in range(0,len(elements)):
			if(elements[ii].type_ == eterna_utils.RNAELEMENT_STACK):
				indices = elements[ii].indices_
				last_pair = sequence[indices[0]] + sequence[indices[1]]

				for jj in range(2,len(indices),2):
					current_pair = sequence[indices[jj]] + sequence[indices[jj+1]]
					if(last_pair == current_pair):
						score -= 1

					last_pair = current_pair

		return score
