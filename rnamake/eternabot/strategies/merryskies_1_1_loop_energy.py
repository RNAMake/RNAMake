import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "1-1 Loop Energy V1"
		self.author_ = "merryskies"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_beta_test_experiment_title_1_1_loop_energy_v1"
		self.default_params_ = [-0.39, 5 ]
		self.code_length_ = 20
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False

	def score(self, design, params):
		elements = design['secstruct_elements']
		sequence = design['sequence']
		pairmap = design['pairmap']

		score = 80
		count = 0
		for ii in range(0,len(elements)):
			elem = elements[ii]
			if(elem.type_ != eterna_utils.RNAELEMENT_LOOP):
				continue

			loop_groups = elem.get_loop_groups()
			closing_pairs = elem.get_loop_closing_pairs(sequence,pairmap)

			if(len(loop_groups) == 2 and len(closing_pairs) == 2):
				if(len(loop_groups[0]) == 1 and len(loop_groups[1]) == 1):
					count += 1
					if(elements[ii].score_ < params[0]):
						score += params[1]


		if(count > 0):
			return score
		else :
			return eterna_utils.UNSCORABLE
