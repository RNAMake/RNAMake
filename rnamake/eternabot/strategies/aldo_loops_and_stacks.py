from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):

		strategy_template.Strategy.__init__(self)

		self.title_ = "Loops & Stacks"
		self.author_ = "aldo"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_loops_stacks"
		self.default_params_ = [2,1,1,2,1,1,0.5, 0.55, 100]
		self.code_length_ = 40
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = True
		self.martin_weight_ = 4.361964
		self.satisfying_point_ = 80

	def score(self, design, params):
		score = 0
		sequence = design['sequence']
		secstruct = design['secstruct']
		elements = design['secstruct_elements']
		total_pairs = design['gc'] + design['gu'] + design['ua']

		for ii in range(0,len(elements)):
			elem = elements[ii]
			#print("AAAA" + str(len(elements)))
			if(elem.type_ == RNAELEMENT_STACK) :
				stack_len = elem.get_stack_length()
				#print("AAAA" + str(stack_len))
				for jj in range(0,stack_len):
					pair = elem.get_pair_from_stack(jj,sequence)
					if(pair == "GC" or pair == "CG"):
						if((jj == 0 or jj == stack_len-1)):
							score += params[0]
						else:
							score += params[1]
					elif(pair == "UA" or pair == "AU"):
						if((jj == 0 or jj == stack_len-1)):
							score += params[2]
						else:
							score += params[3]
					else:
						score += params[4]

		for ii in range(0,len(sequence)):
			if(secstruct[ii] == "." and sequence[ii] == "A"):
				score += params[5]
			elif(secstruct[ii] == "." and sequence[ii] == "G"):
				score += params[6]

		modifier = 1.0 - abs(float(design['gc'])/total_pairs - params[7] )

		#return score
		return params[8] * modifier * (float(score)) / len(sequence)

