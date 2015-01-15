import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):

	def __init__(self):
		strategy_template.Strategy.__init__(self)
		self.title_ = "Berex Test"
		self.author_ = "Berex NZ"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_berex_test"
		self.default_params_ = [0.22, 100, 0.13, 100, 0.2, 100, -60, -30, 1, 97, 107, 1.0]
		self.code_length_ = 30
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = True
		self.martin_weight_ = 10.015664
		self.satisfying_point_ = 70

	def score(self, design, params):
		sequence = design['sequence']
		seqlen = len(sequence)
		seq_range = range(0,len(sequence))
		fe = design['fe']
		mp = design['meltpoint']

		g_count = 0
		c_count = 0
		a_count = 0
		u_count = 0

		for ii in seq_range:
			if(sequence[ii] == "G"):
				g_count += 1
			elif(sequence[ii] == "C"):
				c_count += 1
			elif(sequence[ii] == "A"):
				a_count += 1
			else:
				u_count += 1

		score = 100
		score -= abs(float(g_count)/seqlen - params[0]) * params[1]
		score -= abs(float(u_count)/seqlen - params[2]) * params[3]
		score -= abs(float(c_count)/seqlen - params[4]) * params[5]

		if(fe < params[6]):
			score -= abs(fe - params[6]) * params[8]
		elif(fe > params[7]):
			score -= abs(fe - params[7]) * params[8]

		if(mp < params[9]):
			score -= abs(mp -params[9]) * params[11]
		elif(mp > params[10]):
			score -= abs(mp -params[10]) * params[11]


		return score

