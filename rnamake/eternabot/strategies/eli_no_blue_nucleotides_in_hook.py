import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "No Blue Nucleotides in Hook Area"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_no_blue_nucleotides_in_hook_area"
		self.default_params_ = [-2] # score adjustment to give per blue nucleotide in hook area
		self.code_length_ = 10
		self.publishable_ = True

	def score(self, design, params):
		pairmap = design['pairmap']
		sequence = design['sequence']

		first_pair_index = self.find_nonnegative_value(pairmap)
		last_pair_index = pairmap[first_pair_index] # Last pair index = the pair of the first element in a stack

		score = 100 + (sequence[:first_pair_index+1].count('U') + sequence[last_pair_index:].count('U')) * params[0]
		return score



	def find_nonnegative_value(self, array):
		index = 0
		while(array[index] == -1): index += 1 # Don't need to check for index out of bounds because there are no RNA designs with no stacks?
		return index
