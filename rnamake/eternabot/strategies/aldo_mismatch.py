from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	SHIFT_LIMIT = 3
	#AU GU CG
	PAIR_TYPES = 3
	def __init__(self):
		strategy_template.Strategy.__init__(self)
		self.title_ = "Mismatch"
		self.author_ = "aldo"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_mismatch"
		#ratio of gc pairs in junctions
		self.default_params_ = [40,40,40]
		self.code_length_ = 30
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False
		self.martin_weight_ = -2.793611
		self.satisfying_point_ = -30

	def score(self, design, params):
		sequence = design['sequence']
		pairmap = design['pairmap']
		length = len(pairmap)
		shift = 0
		shift_range = range(0,min(self.SHIFT_LIMIT, length));
		pair_cnt = [0] * self.PAIR_TYPES
		score = 100
		for shift in shift_range:
			i_range = range(0,length-shift)
			for ii in i_range:
				base_a = sequence[ii]
				base_b = sequence[length - ii - shift-1]
				pair_num = self.which_pair(base_a,base_b);
				if(pair_num >= 0):
					pair_cnt[pair_num] += 1
		for ii in range(self.PAIR_TYPES):
			#normalize by sequence length
			score-=(params[ii] * pair_cnt[ii])/length
		return score

	def which_pair(self,base_a, base_b):
		if(base_a>base_b):
			tmp=base_a
			base_a=base_b
			base_b=tmp
		if (base_a == 'A' and base_b=='U'):
			return 0
		elif(base_a == 'G' and base_b == 'U'):
			return 1
		elif(base_a == 'C' and base_b == 'G'):
			return 2
		return -1;
