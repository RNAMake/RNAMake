#by Chengfu chengfuc@andrew.cmu.edu date:7.24
import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template
import re

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)
		self.title_ = "Sameturning GCPair Strategy"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_double_sameturning_gc_pairs"
		self.default_params_ = [1,2]
		self.code_length_ = 30
		self.publishable_ = True

	def score(self, design, params):

		elements = design['secstruct_elements']
		sequence = design['sequence']

		score = 0

		first_stack = True
		for ii in range(0,len(elements)):
			if(elements[ii].type_ == eterna_utils.RNAELEMENT_STACK):
				last_pair = ""
				sameturning = 0
				previous_sameturned = False
				consecutive_turns = 0
				for jj in range(0,len(elements[ii].indices_),2):
					pair = sequence[elements[ii].indices_[jj]] + sequence[elements[ii].indices_[jj+1]]
					pair = pair.upper()

					if(pair == "GC" or pair == "CG"):
						if(pair == last_pair):
							sameturning += 1

							if(first_stack):
								if(sameturning > 2):
									score -= params[0]
								elif(previous_same_turned):
									score -= params[0]
							else:
								if(len(elements[ii].indices_) >= 12):
									if(sameturning > 1):
										score -= params[1]
								else:
									score -= params[1]

							previous_same_turned = True
						else:
							previous_same_turned = False
					else:
						previous_same_turned = False

					last_pair = pair
				first_stack = False
		return score









