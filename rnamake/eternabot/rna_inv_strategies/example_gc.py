import random
import inv_utils
class Strategy:
	def __init__(self):
		self.title_ = "[Example] Randomly fill gc pairs"
		self.author_ = "Example"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_example_60_of_pairs_must_be_gc_pairs-1erc6"
		self.code_length_ = 1
		self.publishable_ = True
		
	def solve(self, design):
		bases=inv_utils.BASES
		random.seed()
		elements = design['secstruct_elements']
		target_pair = design['secstruct']
		pair_map = design['pairmap']
		length = len(target_pair)
		sequence = ['A']*length
		natural_pair = ""
		while(natural_pair != target_pair):
			for elem in elements:
				inv_utils.fill_gc(elem, pair_map , sequence, random)
			ret = inv_utils.fold(sequence)
			natural_pair = ret[0]
		return sequence