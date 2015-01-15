from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template


class Strategy(strategy_template.Strategy):
	def __init__(self):
		self.title_ = "NUPACK test"
		self.author_ = "NUPACK"
		self.url_ = "http://eterna.cmu.edu/htmls/player.html?profid=24553"
		self.default_params_ = []
		self.code_length_ = 1
		self.publishable_ = True


	def score(self, design, params):
		return 100.0 - float(design['nupack_ed']) * 20;

