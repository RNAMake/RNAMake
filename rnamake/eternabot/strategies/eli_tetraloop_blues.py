from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Tetraloop blues"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_tetraloop_blues"
		#ratio of gc pairs in junctions
		self.default_params_ = [4]
		self.code_length_ = 30
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False


	def score(self, design, params):
		sequence = design['sequence']
		elements = design['secstruct_elements']
		pairmap = design['pairmap']

		penalty= 0
		count = 0
		for ii in range(0,len(elements)):
			if(elements[ii].type_ == RNAELEMENT_LOOP
				and len(elements[ii].children_)==0
				and elements[ii].parent_
				and len(elements[ii].indices_)==4):
				count += 1
				indices=elements[ii].indices_
				u_count=0
				for jj in range(0,len(indices)):
					if(sequence[indices[jj]]=="U"):
						u_count += 1

				#pair at the closing stack
				parent = elements[ii].parent_
				npi = len(parent.indices_)
				for kk in range(1,3):
					if(sequence[parent.indices_[npi-kk]]=="U"):
						u_count += 1

				if(u_count>2):
					penalty += u_count-2
		if(count == 0):
			return UNSCORABLE
		score = 100 -penalty*params[0];
		return score
