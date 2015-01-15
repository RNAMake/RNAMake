from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "GC-pairs in junctions"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_gc_pairs_in_junctions"
		#ratio of gc pairs in junctions
		self.default_params_ = [1]
		self.code_length_ = 30
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False
		self.satisfying_point_ = 0

	def score(self, design, params):
		sequence = design['sequence']
		pairmap = design['pairmap']
		length = len(pairmap)
		i_range = range(0,length)
		junct_cnt = 0
		#number of gc pairs in juctions
		junct_gc = 0
		score = 100
		for ii in i_range:
			idx=pairmap[ii]
			#do not double count
			if (idx>ii):
				#check neighbors
				if (not self.is_neibor_paired(ii,pairmap,length) )\
					or (not self.is_neibor_paired(idx,pairmap,length)):
					junct_cnt += 1
					if(self.is_gc(sequence[ii],sequence[idx])):
						junct_gc += 1

                score += - params[0] * (junct_cnt - junct_gc)

		return score

	def is_neibor_paired(self,idx, pairmap,length):
		if( (idx>0) and (pairmap[idx-1]<0) ) \
			or ( (idx<length-1) and (pairmap[idx+1]<0) ):
			return False
		return True

	def is_gc(self, base_a, base_b):
		if(base_a=="G" and base_b=="C") \
			or(base_a=="C" and base_b=="G"):
			return True
		return False
