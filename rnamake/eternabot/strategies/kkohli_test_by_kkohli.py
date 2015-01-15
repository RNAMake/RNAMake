import rnamake.eternabot.eterna_utils as eterna_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Test by kkohli"
		self.author_ = "kkohli"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_test_by_kkohli"
		self.default_params_ = [1.5, 10, 10, 1, 10, 10]
		self.code_length_ = 30
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = True

	def score(self, design, params):

		elements = design['secstruct_elements']
		meltpoint = design['meltpoint']
		sequence = design['sequence']
		pairmap = design['pairmap']
		total_pairs = design['gc'] + design['gu'] + design['ua']
		fe = design['fe']
		score = 100

		g_count = 0
		for ii in range(0,len(sequence)):
			if(sequence[ii] == "G"):
				g_count += 1
		g_rate = float(g_count) / len(sequence)

		closing_pairs_count = 0
		closing_gc_count = 0
		for ii in range(0,len(elements)):
			if(elements[ii].type_ == eterna_utils.RNAELEMENT_LOOP):
				if(elements[ii].score_ > params[0]):
					score -= params[1]

				closing_pairs = elements[ii].get_loop_closing_pairs(sequence,pairmap)
				for jj in range(0,len(closing_pairs)):
					closing_pairs_count += 1
					if(closing_pairs[jj] == "GC" or closing_pairs[jj] == "CG"):
						closing_gc_count += 1
		closing_gc_rate = float(closing_gc_count) / closing_pairs_count


		if(meltpoint < 97 or meltpoint > 107):
			score -= params[2]

		if(fe > -30):
			score -= (fe - (-30)) * params[3]

		if(g_rate < 0.25 or g_rate > 0.35):
			score -= params[4]

		if(closing_gc_rate < 0.6):
			score -= (0.6 - closing_gc_rate) * params[5]

		return score

