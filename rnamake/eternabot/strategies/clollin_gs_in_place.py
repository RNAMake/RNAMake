import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)
		self.title_ = "G's in place of the last A's on the right hand side of any end loop"
		self.author_ = "clollin"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_ensemble_strategy_or_beta_test_experiment_title_gs_in_place_of_the_last_as_on_the_right"
		self.default_params_ = []
		self.code_length_ = 20
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False

	def score(self, design, params):

		sequence = design['sequence']
		elements = design['secstruct_elements']
		count = 0
		score = 80

		for ii in range(0,len(elements)):
			if (elements[ii].type_ == eterna_utils.RNAELEMENT_LOOP):
				loop_groups = elements[ii].get_loop_groups()

				# For hairpin with loop longer than 3
				if(len(loop_groups) == 1 and len(loop_groups[0]) > 3):
					count +=1
					if(sequence[loop_groups[0][0]] == "G"):
						score += 1


		if(count < 1):
			return eterna_utils.UNSCORABLE

		score = (score)

		return score
