import  rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):

	def __init__(self):
		strategy_template.Strategy.__init__(self)
		self.title_ = "Simplified Berex Test"
		self.author_ = "jeehyung"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_simplified_berex_test"
		self.default_params_ = [1,1, -60, -30]
		self.code_length_ = 30
		self.publishable_ = True

	def score(self, design, params):
		sequence = design['sequence']
		seqlen = len(sequence)
		seq_range = range(0,len(sequence))
		fe = design['fe']
		mp = design['meltpoint']

		score = 100

		if(fe < -60):
			score -= abs(-60.0 - fe) * params[0]
		elif(fe > -30):
			score -= abs(-30.0 - fe) * params[0]

		if(mp < 97):
			score -= abs(97.0 - mp) * params[1]
		elif(mp > 107):
			score -= abs(107.0 - mp) * params[1]

		return score

	def patch(self, design, params):
		#Change GCs to AUs from the stack with the least GC percentage
		return design
