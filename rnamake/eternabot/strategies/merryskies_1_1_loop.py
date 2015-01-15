from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "The 1-1 Loop V1"
		self.author_ = "merryskies"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_beta_test_experiment_title_the_1_1_loop_v1"
		self.default_params_ = [5,3,1]
		self.code_length_ = 20
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False

	def score(self, design, params):
		elements = design['secstruct_elements']
		sequence = design['sequence']
		pairmap = design['pairmap']

		score = 80

		for ii in range(0,len(elements)):
			elem = elements[ii]
			if(elem.type_ != RNAELEMENT_LOOP):
				continue

			loop_groups = elem.get_loop_groups()
			closing_pairs = elem.get_loop_closing_pairs(sequence,pairmap)

			if(len(loop_groups) == 2 and len(closing_pairs) == 2): #case 2,3,4
				if(len(loop_groups[0]) == 1 and len(loop_groups[1]) == 1): #case 2
					if(sequence[loop_groups[0][0]] == "G" and sequence[loop_groups[1][0]] == "G"):
						score += params[0]
					elif(sequence[loop_groups[0][0]] == "G" and sequence[loop_groups[1][0]] == "A"):
						score += params[1]
					elif(sequence[loop_groups[0][0]] == "A" and sequence[loop_groups[1][0]] == "G"):
						score += params[1]
		 			elif(sequence[loop_groups[0][0]] == "A" and sequence[loop_groups[1][0]] == "A"):
						score += params[2]


		return score
