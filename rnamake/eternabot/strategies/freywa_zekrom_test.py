import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Zekrom Test"
		self.author_ = "Freywa"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_zekrom_test"
		self.default_params_ = [0.5, 100, 2]
		self.code_length_ = 10
		self.publishable_ = True


	def score(self, design, params):
		gcs = design['gc']
		gus = design['gu']
		aus = design['ua']
		fe = design['fe']
		npairs = gcs + gus + aus

		score = 0

		gcp = float(gcs)/npairs
		score -= abs(gcp - params[0]) * params[1]

		score -= gus * params[2]

		return score

