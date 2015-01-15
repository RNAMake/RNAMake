from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)
		self.title_ = "Berex Loop Basic"
		self.author_ = "Berex NZ"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_berex_loop_basic"
		self.default_params_ = [1,3,2,3,3,3,1,2]
		self.code_length_ = 70
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
			indices =elem.indices_

			if(len(loop_groups) == 1 and len(closing_pairs) == 2): #case 1
				if(closing_pairs[0] == "GC" or closing_pairs[0] == "CG" or closing_pairs[1] == "GC" or closing_pairs[1] == "CG"):
					score += params[0]
			elif(len(loop_groups) == 2 and len(closing_pairs) == 2): #case 2,3,4
				if(len(loop_groups[0]) == 1 and len(loop_groups[1]) == 1): #case 2
					if(sequence[loop_groups[0][0]] == "G" and sequence[loop_groups[1][0]] == "G"):
						score += params[1]
				elif(len(loop_groups[0]) + len(loop_groups[1]) == 3): #case 3
					if(len(loop_groups[0]) == 1):
						bottom = 0
						top = 1
					else:
						bottom = 1
						top = 0

					if(sequence[loop_groups[bottom][0]] == "G"):
						score += params[2]
				elif(len(loop_groups[0]) + len(loop_groups[1]) == 4): #case 4
					if(len(loop_groups[0]) == 2 and len(loop_groups[1]) == 2):
						if((sequence[loop_groups[0][0]] == "G" or sequence[loop_groups[0][1]] == "G") and (sequence[loop_groups[1][0]] == "G" or sequence[loop_groups[1][1]] == "G") ):
							score += params[3]
					else:
						if(len(loop_groups[0]) == 1):
							bottom = 0
							top = 1
						else:
							bottom = 1
							top = 0

						top_g_count = 0
						if(sequence[loop_groups[top][0]] == "G"):
							top_g_count += 1

						if(sequence[loop_groups[top][1]] == "G"):
							top_g_count += 1

						if(sequence[loop_groups[top][2]] == "G"):
							top_g_count += 1

						if(sequence[loop_groups[bottom][0]] == "G"):
							score += params[4]

						if(top_g_count == 2):
							score += params[4]

				else : #case 5
					if(len(loop_groups) == 1 and len(loop_groups[0]) == 4 and len(closing_pairs) == 1):
						if(sequence[loop_groups[0][0]] == "G"):
							score += params[5]
					else:
						g_count = 0
						for jj in range(0,len(indices)):
							if(sequence[indices[jj]] == "G"):
								g_count+= 1

						if(g_count == 1):
							score += params[6]
						elif(g_count >=2 ):
							score += params[7]

		return score

