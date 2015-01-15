import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Tests by Region - boundaries"
		self.author_ = "Quasispecies"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_tests_by_region"
		#ratio of gc pairs in junctions
		self.default_params_ = [1]
		self.code_length_ = 10
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False

	def score(self, design, params):

		score = 100
		sequence = design['sequence']
		elements = design['secstruct_elements']

		for ii in range(0,len(elements)):
			if(elements[ii].type_ == eterna_utils.RNAELEMENT_STACK):
				if(len(elements[ii].indices_) >= 4):
					if(elements[ii].get_pair_from_stack(0,sequence) == elements[ii].get_pair_from_stack(1,sequence)):
						score -= self.default_params_[0]
				if(len(elements[ii].indices_) >= 6):
					last_pair = len(elements[ii].indices_) / 2 -1
					if(elements[ii].get_pair_from_stack(last_pair,sequence) == elements[ii].get_pair_from_stack(last_pair-1,sequence)):
						score -= self.default_params_[0]
		return score
