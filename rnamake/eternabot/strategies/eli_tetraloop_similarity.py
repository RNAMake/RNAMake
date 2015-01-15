from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Tetraloop similarity"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_tetraloop_similarity"
		#ratio of gc pairs in junctions
		self.default_params_ = [4,4, 0.3]
		self.code_length_ = 40
		self.publishable_ = True
		#energies within this difference are treated as similar

	def score(self, design, params):
		sequence = design['sequence']
		elements = design['secstruct_elements']
		pairmap = design['pairmap']

		penalty= 0
		count = 0
		tetraloop_idx = []
		energies = []
		for ii in range(0,len(elements)):
			if(elements[ii].type_ ==RNAELEMENT_LOOP
				and len(elements[ii].children_)==0
				and elements[ii].parent_
				and len(elements[ii].indices_)==4):
				count += 1
				tetraloop_idx.append(ii)
				energies.append(elements[ii].score_)
		score = 0;
		if(count == 4):
			energies.sort()
			#all energies are similar
			if(energies[3]-energies[0]< params[2]):
				score += 5
			#two pairs
			elif(energies[3]-energies[2] < params[2]
				and energies[1] - energies[0] < params[2]):
				score += 5
			#three similar loops
			elif(energies[3]-energies[1]< params[2]
				or energies[2]-energies[0]< params[2]):
				score += 3
			elif(energies[3]-energies[2] < params[2] or
				energies[2]-energies[1] < params[2] or
				energies[1]-energies[0] < params[2]):
				score+=2

			return score*20
		elif(count==3):
			energies.sort()
			if(energies[2]-energies[0] < params[2]):
				score+=5
			elif(energies[2]-energies[1] < params[2]
				or energies[1]-energies[0] < params[2]):
				score+=3
			return score*20
		else:
			return UNSCORABLE
