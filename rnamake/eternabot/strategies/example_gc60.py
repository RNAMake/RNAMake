import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):

		strategy_template.Strategy.__init__(self)

		self.title_ = "[Example] 60% of pairs should be GC pairs"
		self.author_ = "Example"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_example_60_of_pairs_must_be_gc_pairs-1erc6"
		self.default_params_ = [0.6]
		self.code_length_ = 1
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive = False


	def score(self, design, params):
		return 100-(abs(params[0] - float(design['gc']) / (design['gu'] + design['gc'] + design['ua']))) * 100
