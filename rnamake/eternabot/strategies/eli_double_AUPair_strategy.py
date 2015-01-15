#by Chengfu chengfuc@andrew.cmu.edu date:6.19
from rnamake.eternabot.eterna_utils import *
import re
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Double AU-pair strategy"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_double_au_pair_strategy"
		self.default_params_ = [-1]
		self.code_length_ = 190
		self.publishable_ = False


	#return the index of neck pair and the neck pair
	def get_neck_pairs(self,loop_groups,sequence,pairmap):

		pairs = []
		pair_indices = []


		for ii in range(0,len(loop_groups)):
			start_i = loop_groups[ii][0]
			end_i = loop_groups[ii][len(loop_groups[ii]) - 1]

			if(start_i - 1 >= 0):
				if(pairmap[start_i-1] <0):
					print "ERROR something's wrong with RNAELEMENT 1"
					sys.exit(0)

				if(pair_indices.count(start_i-1) == 0):
					pair = sequence[min(start_i-1,pairmap[start_i-1])] + sequence[max(start_i-1,pairmap[start_i-1])]
                                        #getting the index and pair
					pairs.append(str(min(start_i-1,pairmap[start_i-1]))+pair.upper())
					pair_indices.append(start_i-1)
					pair_indices.append(pairmap[start_i-1])

			if(end_i + 1 < len(sequence)):
				if(pairmap[end_i+1] <0):
					print "ERROR something's wrong with RNAELEMENT 2"
					sys.exit(0)

				if(pair_indices.count(end_i+1) == 0):
					pair = sequence[min(end_i+1,pairmap[end_i+1])] + sequence[max(end_i+1,pairmap[end_i+1])]
					pairs.append(str(min(end_i+1,pairmap[end_i+1]))+pair.upper())
					pair_indices.append(end_i+1)
					pair_indices.append(pairmap[end_i+1])

		return pairs


        def score(self, design, params):

		elements = design['secstruct_elements']
		sequence = design['sequence']
		pairmap = design['pairmap']
		neckArea = []
		pairStackGrps = []
                adjPairCounter=0
		score = 100
		#================================================================================================Get necks========================
		for ii in range(0,len(elements)):                                                                                               #=
			elem = elements[ii]
			if(elem.type_ == RNAELEMENT_LOOP):
                                loop_groups = elem.get_loop_groups()
                                tmp=[]
                                tmp = self.get_neck_pairs(loop_groups,sequence,pairmap)
                                for jj in range(0, len(tmp)):
                                        neckArea.append(tmp[jj])

                #print(neckArea)#=
                #=================================================================================================================================


                #Getting all the adjacent au pairs================================================Adjacent AuPairs================================
                for ii in range(0,len(elements)):                                                                                               #=
                        elem = elements[ii]                                                                                                     #=
                        if(elem.type_ == RNAELEMENT_STACK):
                                                                                                                                                #=
                                lengthOfStack = elem.get_stack_length()
                                pairStackRaw = []                                                                                               #=
                                indexOfPair = []
                                rawIndex = []
                                #counter=0
                                for jj in range(0,lengthOfStack):
                                        rawIndex.append(-1)
                                for jj in range(0,lengthOfStack):
                                        pair = elem.get_pair_from_stack(jj,sequence)                                                            #=

                                        if(pair == "UA" or pair == "AU"):                                                                       #=
                                                indexOfPair.append(jj)
                                                #rawIndex[jj] =  counter
                                                #counter = counter +1
                                                pairContent = (sequence[elem.indices_[jj * 2]] + sequence[elem.indices_[jj * 2 + 1]]).upper()
                                                pairContent = str(min(elem.indices_[jj * 2],elem.indices_[jj * 2 + 1]))+pairContent
                                                pairStackRaw.append(pairContent)
                                counter=0
                                for jj in range(1,len(indexOfPair)):
                                        if(indexOfPair[jj]-indexOfPair[jj-1]==1):                                                               #=
                                                stack=[]                                                                                        #=
                                                stack.append(pairStackRaw[jj-1])
                                                stack.append(pairStackRaw[jj])
                                                pairStackGrps.append(stack)
                                #print(pairStackGrps)
                                #if(counter>0):                                                                                                  #=
                                #        score += (counter-1)*params[0];
                ############======================================================================================================================
                ############=========================================================================Now Comparing Neck and Pair==================
                ############======================================================================================================================

                #caseA=========================== waive
                #caseB=========================== penalty

                #if there is 2 double au pair(4 pairs and twisted) in the neck, waive the penalty-----------------------------------------caseA
                checker = 0
                tmp=[]
                for kk in range(1,len(pairStackGrps)):
                        n1=re.findall(r"\d",pairStackGrps[kk-1][1])
                        n2=re.findall(r"\d",pairStackGrps[kk][0])
                        #if 4 in row
                        if(n1!=n2):
                                continue
                        #if twisted
                        if(re.findall(r"\d",pairStackGrps[kk-1][0])!=re.findall(r"\d",pairStackGrps[kk-1][1]) and\
                           re.findall(r"\d",pairStackGrps[kk][0]) != re.findall(r"\d",pairStackGrps[kk][1]) and\
                           re.findall(r"\d",pairStackGrps[kk-1][1]) != re.findall(r"\d",pairStackGrps[kk][0])):
                                print("this never happen")
                                for ss in range(0,len(neckArea)):
                                    if(neckArea[ss]==pairStackGrps[kk-1][0] or neckArea[ss]==pairStackGrps[kk][1]):
                                            return 100#============================================================without penalty but never happened
                                pairStackGrps.pop(kk-1)
                                pairStackGrps.pop(kk)

                #--------------------------------------------------------------------------------------------------------------------End of CaseA
                ############======================================================================================================================
                #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                #au pairs in the neck are not counting&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&caseB1
                #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


                pairStackTuningSame=[]

                for kk in range(0,len(pairStackGrps)):
                        for mm in range(1, len(pairStackGrps[kk])):
                                a=re.findall(r"\D",pairStackGrps[kk][mm-1])
                                b=re.findall(r"\D",pairStackGrps[kk][mm])
                                if(a==b):
                                        stack=[]
                                        stack.append(pairStackGrps[kk][mm-1])
                                        stack.append(pairStackGrps[kk][mm])
                                        pairStackTuningSame.append(stack)

                adjPairCounter = len(pairStackTuningSame)
                #print(pairStackTuningSame,adjPairCounter)

                for kk in range(0,len(pairStackTuningSame)):
                        numOfNeck=0
                        for mm in range(0, len(pairStackTuningSame[kk])):
                                for yy in range(0, len(neckArea)):
                                        if(pairStackTuningSame[kk][mm]==neckArea[yy]):#&&&&&&&&&&&&&&&&&&Satisfy caseB1, del grp
                                                numOfNeck=numOfNeck+1
                        if(numOfNeck>0):
                                #print("@@@@@Neck",adjPairCounter)
                                adjPairCounter=adjPairCounter-1
                #----&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&End of CaseB1

                ############======================================================================================================================
                if(adjPairCounter>0):
                        score += (adjPairCounter-1)*params[0]#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%More than One Double Pairs
                return score
















