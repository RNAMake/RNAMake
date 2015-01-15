import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Xmas, cubs and illusions"
		self.author_ = "Edward Lane"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_xmas_cubs_and_illusions"
		self.default_params_ = []
		self.code_length_ = 10
		self.publishable_ = True


	def score(self, design, params):
		gcs = design['gc']
		gus = design['gu']
		aus = design['ua']

		score = 100

		if(gus + aus == 0):
			score -= 100
		elif(gcs + aus == 0):
			score -= 100
		elif(gcs + gus == 0):
			score -= 100

		return score
