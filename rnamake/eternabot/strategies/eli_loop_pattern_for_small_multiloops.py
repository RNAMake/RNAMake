from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Loop pattern for small multiloops"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_loop_pattern_for_small_loops"
		self.default_params_ = [100, -1, -2] # [Score for Good GC multiloop, Amount to modify good score if a GC pair is wrong direction, Score for Non-GC pair in multiloop]
		self.code_length_ = 65
		self.publishable_ = True
		self.skippable_indices = []
		self.good_patterns = ["CGGCGCCG", "GCCGCGGC", "CGCGCGCG"] # These are the patterns that Eli specified as ideal


	def score(self, design, params):
                pairmap = design['pairmap']

                score = params[0]
                gc_multiloop_found = False
                self.skippable_indices = []

                for ii in range(len(pairmap)):
                        if ii in self.skippable_indices: continue
                        gc_multiloop_string = self.get_gc_multiloop_members(ii, design)
                        if gc_multiloop_string != "":
                                if not gc_multiloop_found: gc_multiloop_found = True

                                # The following lines loop through the multiloop and record which indices we've already explored so we can skip them later.
                                self.skippable_indices.append(ii) # Current index
                                self.skippable_indices.append(pairmap[ii]) # Current index's pair
                                ''' # These are commented out to make sure we don't exclude any funky connected multi-loops that might be designed in the future.
                                self.skippable_indices.append(pairmap[ii]-1) # Pair Index-1
                                self.skippable_indices.append(pairmap[pairmap[ii]-1]) # Pair
                                self.skippable_indices.append(pairmap[pairmap[ii]-1]-1) # Pair Index-1
                                self.skippable_indices.append(pairmap[pairmap[pairmap[ii]-1]-1]) # Pair
                                self.skippable_indices.append(pairmap[pairmap[pairmap[ii]-1]-1]-1)# Pair Index-1
                                self.skippable_indices.append(pairmap[pairmap[pairmap[pairmap[ii]-1]-1]-1]) # Pair
                                '''

                                # Modify score for each error in gc loop.
                                if gc_multiloop_string in self.good_patterns: continue
                                else:
                                        errors = self.compare_pairs(self.good_patterns[0], gc_multiloop_string, params)
                                        errors = max(errors, self.compare_pairs(self.good_patterns[1], gc_multiloop_string, params))
                                        errors = max(errors, self.compare_pairs(self.good_patterns[2], gc_multiloop_string, params))
                                        score += errors
                if not gc_multiloop_found: return UNSCORABLE
		return score # Returned if there are no multiloops with GC pairs only

        # Returns the number of differences between two strings of 4 pairs.
        # First string is the "good" string for which we will compare the 2nd to for errors.
	def compare_pairs(self, string1, string2, params):
                errors = 0
                if string1[0:2] != string2[0:2]:
                        if (self.has_non_gc(string2[0:2])): errors += params[2]
                        else: errors += params[1]
                if string1[2:4] != string2[2:4]:
                        if (self.has_non_gc(string2[2:4])): errors += params[2]
                        else: errors += params[1]
                if string1[4:6] != string2[4:6]:
                        if (self.has_non_gc(string2[4:6])): errors += params[2]
                        else: errors += params[1]
                if string1[6:8] != string2[6:8]:
                        if (self.has_non_gc(string2[6:8])): errors += params[2]
                        else: errors += params[1]
                return errors

	# This function returns a string of the multiloop members
	# returns an empty string if there is an error (ie. not a multiloop, already looped over a pair)
	# As soon as we find a GC multiloop, we know we're viewing it from the neck (lowest index). We'll skip it otherwise with self.skippable_indices.
        # Strings go in this order: NECK BOTTOM LEFT TOP (clockwise), with first element being LEFT and second being RIGHT
	def get_gc_multiloop_members(self, index, design):
                returnstring = ""
                sequence = design['sequence']
                pairmap = design['pairmap']
                paired_index = pairmap[index]
                returnstring += sequence[index]
                if (paired_index != -1):
                        returnstring += sequence[paired_index]
                        index2 = paired_index-1
                        if index2 in self.skippable_indices: return ""
                        returnstring += sequence[index2]
                        paired_index2 = pairmap[index2]
                        if (paired_index2 != -1):
                                returnstring += sequence[paired_index2]
                                index3 = paired_index2-1
                                if (index3 in self.skippable_indices) or (index3 == index): return "" # If we've done this loop we don't need to do it again, and make sure we're not just looping over 2 stacks
                                returnstring += sequence[index3]
                                paired_index3 = pairmap[index3]
                                if (paired_index3 != -1):
                                        returnstring += sequence[paired_index3]
                                        index4 = paired_index3-1
                                        if (index4 in self.skippable_indices) or (index4 == index): return "" # If we've done this loop we don't need to do it again, and make sure we're not looping over just 3 stacks
                                        returnstring += sequence[index4]
                                        paired_index4 = pairmap[index4]
                                        if (paired_index4 != -1):
                                                returnstring += sequence[paired_index4]
                                                return returnstring
                return ""

        # Checks if there are any non-G or non-C elements in a string.
        def has_non_gc(self, string):
                if (string.find("A") != -1): return True
                if (string.find("U") != -1): return True
                return False
