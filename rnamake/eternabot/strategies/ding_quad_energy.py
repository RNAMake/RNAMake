#by Chengfu chengfuc@andrew.cmu.edu date:6.27
from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)
		self.title_ = "Ding quad energy"
		self.author_ = "Ding"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/tell_us_about_your_eterna_lab_algorithms"
                #0: energies over -0.9
		#1: penalty -1
		#------------------
		#2: energies -2.1
		#3: energies -2.5
		#4: penalty -1
		#------------------
		#5: penalty for 3 AU -1
		#6: 10 defining for langer stem
		#------------------
		#7: 2 stems closed with AU
		#8: penalty -1
		self.default_params_ = [-0.9,-1,-2.10,-2.50,-1,-1,10,2,-1]
		self.code_length_ = 190
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = False

	def score(self, design, params):

		elements = design['secstruct_elements']
		sequence = design['sequence']
		pairmap = design['pairmap']


		score = 100
		GCEndingStems=0
		AUEndingStems=0
		TotalStems=0
                for ii in range(0,len(elements)):
			elem = elements[ii]

			if(elem.type_ == RNAELEMENT_STACK):
                                TotalStems=TotalStems+1

				for ii in range(0,len(elem.quad_scores_)):
                                        #params[0] = -0.9
                                        if(elem.quad_scores_[ii]/100.0>params[0]):
                                                score+=params[1]
                                lengthOfStack = elem.get_stack_length()

                                for ii in range(0, lengthOfStack-1):
                                        pair = elem.get_pair_from_stack(ii,sequence)
                                        nextPair= elem.get_pair_from_stack(ii+1,sequence)

                                        if(pair == "GU" or pair == "UG" or nextPair == "GU" or nextPair == "UG"):
                                           if(elem.quad_scores_[ii]/100.0 > params[2]):
                                                   score += params[4]

                                if(lengthOfStack<params[6]):
                                        for ii in range(0, lengthOfStack-2):
                                                pair0=elem.get_pair_from_stack(ii,sequence)
                                                pair1=elem.get_pair_from_stack(ii+1,sequence)
                                                pair2=elem.get_pair_from_stack(ii+2,sequence)

                                                if((pair0=='AU' or pair0=='UA') and\
                                                   (pair1=='AU' or pair1=='UA') and\
                                                   (pair2=='AU' or pair2=='UA')):
                                                        score+=params[5]

                                pairStart=elem.get_pair_from_stack(0,sequence)
                                pairEnd=elem.get_pair_from_stack(lengthOfStack-1,sequence)

                                if((pairStart=='GC' or pairStart=='CG') and\
                                   (pairEnd=='GC' or pairEnd=='CG')):
                                        GCEndingStems = GCEndingStems+1
                                if((pairStart=='AU' or pairStart=='UA') and\
                                   (pairEnd=='AU' or pairEnd=='UA')):
                                        AUEndingStems = AUEndingStems+1
                multiplier=TotalStems-GCEndingStems
                if(multiplier>0 and AUEndingStems-params[7]>0):
                        score+=params[8]*(multiplier- (min(AUEndingStems,params[7])))

                return score
















