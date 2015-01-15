import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Clean Dot Plot"
		self.author_ = "penguian"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_clean_dot_plot"
		self.default_params_ = []
		self.code_length_ = 10
		self.publishable_ = True
		self.is_part_of_ensemble_ = False
		self.denormalized_ = True
		self.comprehensive_ = False

	def score(self, design, params):

		penalty = 0.0

		n = len(design['sequence'])
		npairs = design['gc'] + design['gu'] + design['ua']
		dotplot = design['dotplot']
		pairmap = design['pairmap']

		for ii in range(0,len(dotplot)):
			i_index = dotplot[ii][0]
			j_index = dotplot[ii][1]

			if(pairmap[i_index] != j_index):
				penalty += dotplot[ii][2]

		return 100 -(penalty / npairs)

