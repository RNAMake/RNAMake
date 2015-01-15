#by Chengfu chengfuc@andrew.cmu.edu date:6.11
from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template


class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "JP-Stratmark"
		self.author_ = "JerryP70"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_jp_stratmark"
		#0:Score for Tetraloops and pentaloops and any other larger loops
		#1:the number using for ckecking pair continuity
		#2:the number using for checking bonded pairs
		#3: 45%
		#4: 65%
		#5: 40%
		#6: 50%
		#7: 5%

		self.default_params_ = [10,4,6,0.45,0.65,0.4,0.5,0.05,10]
		self.code_length_ = 150
		self.publishable_ = True
		self.denormalized_ = True
		self.comprehensive_ = True

	#Get the stacks pair nearing the LoopGrps
	#by given the numOfStacks for recursively searching
	#return stacks containing the looping index and stacks informations
	def GetClosingStacks(self,loopGrps,sequence,pairmap,numOfStacks):

		stacks = []

		#for storing the root pair indices
		rootPairIndices = []
                #print('GetClosing')
                #print(loopGrps)
                for ii in range(0, len(loopGrps)):
                        startIndex = loopGrps[ii][0]
                        if(startIndex - 1 >= 0):

                                #Making sure it is bonded pair
                                if(pairmap[startIndex-1] <0):
                                        return []

                                #Add pairs into stacks
                                if(rootPairIndices.count(startIndex-1) == 0):
                        		rootPairIndices.append(startIndex-1)
                                        rootPairIndices.append(pairmap[startIndex-1])
                                        pairs = []
                                        #writing the loop index
                                        pairs.append(str(ii))

                                        for jj in range(1,numOfStacks+1):

                                                if(pairmap[startIndex-jj]<0):
                                                        break

                                                pair = sequence[min(startIndex-jj,pairmap[startIndex-jj])] + sequence[max(startIndex-jj,pairmap[startIndex-jj])]
                                                #writing pair
                                                pairs.append(pair.upper())

                                        #end of 'for' in enumerating numOfStacks
                                        stacks.append(pairs)


                        endIndex = loopGrps[ii][len(loopGrps[ii]) - 1]
                        if(endIndex + 1 < len(sequence)):
                                #Making sure it is bonded pair
                                if(pairmap[endIndex+1] <0):
                                        return []
                                #Add pairs into stacks
                                if(rootPairIndices.count(endIndex+1) == 0):
                        		rootPairIndices.append(endIndex+1)
                                        rootPairIndices.append(pairmap[endIndex+1])
                                        pairs = []
                                        #writing the loop index
                                        pairs.append(str(ii))

                                        for jj in range(1,numOfStacks+1):

                                                if(pairmap[endIndex+jj]<0):
                                                        break

                                                pair = sequence[min(endIndex+jj,pairmap[endIndex+jj])] + sequence[max(endIndex+jj,pairmap[endIndex+jj])]
                                                #writing pair
                                                pairs.append(pair.upper())

                                        #end of 'for' in enumerating numOfStacks
                                        stacks.append(pairs)

                #end of 'for' in loops
                return stacks


        def score(self, design, params):

		elements = design['secstruct_elements']
		sequence = design['sequence']
		pairmap = design['pairmap']

		score = 0
                for ii in range(0,len(elements)):
			elem = elements[ii]
			if(elem.type_ == RNAELEMENT_LOOP):

                                loop_groups = elem.get_loop_groups()
                                closing_pairs = elem.get_loop_closing_pairs(sequence,pairmap)

                                #tatraloop and pentaloops:
                                if(len(loop_groups) == 1 and len(closing_pairs) == 1):
                                        lengthOfGrpIndices = len(loop_groups[0])

                                        if(lengthOfGrpIndices==4 or lengthOfGrpIndices==5):
                                                if(closing_pairs[0] == "GC" or closing_pairs[0] == "CG"):
                                                        score += params[0] #*lengthOfGrpIndices
                                                        #print("Tatraloops and Pentaloops:"+str(params[0]*lengthOfGrpIndices))


                                #more larger loops:
                                totalNum = 0
                                for jj in range(0,len(loop_groups)):
                                        for kk in range(0,len(loop_groups[jj])):
                                                totalNum = totalNum+1

                                #skip if don't has large loop
                                #if(!checker):
                                if(totalNum<5):
                                        continue
                                #Getting 2 adjacent stacks nearby large loop
                                stacks = self.GetClosingStacks(loop_groups,sequence,pairmap,2)
                                #could be error given by GetClosingStacks
                                if(len(stacks)==0):
                                        continue

                                for jj in range(0, len(stacks)):
                                        if (len(stacks[jj])<3):
                                                continue
                                        #stacks[][0] is loop index
                                        if((stacks[jj][1]=="GC" or stacks[jj][1]=="CG") and (stacks[jj][2]=="GC" or stacks[jj][2]=="CG")):
                                                multiplier = len(loop_groups[int(stacks[jj][0])])
                                                score += params[0] #*multiplier
                                                #print("LargerLoops"+str(params[0]*multiplier))

                        if(elem.type_ == RNAELEMENT_STACK):

                                lengthOfStack = elem.get_stack_length()
                                #params[2] is 6
                                #skip when length is less than 6
                                if(lengthOfStack<params[2]):
                                        continue
                                indexOfPair = []
                                for jj in range(0,lengthOfStack):
                                        pair = elem.get_pair_from_stack(jj,sequence)
					if(pair == "GC" or pair == "CG"):
                                                indexOfPair.append(jj);
                                checker = 1
                                for jj in range(1,len(indexOfPair)):
                                        #params[1] is 4
                                        if(indexOfPair[jj]-indexOfPair[jj-1]>=params[1]):
                                                checker = 0
                                                break
                                if(checker):
                                        score += params[0]; #*lengthOfStack
                                        #print("ContinueStack:"+str(params[0]*lengthOfStack))


                #end of for every element

                #params[3]: 45% gc
                #params[4]: 65% gc
                #params[5]: 40% ua
                #params[6]: 50% ua
                #params[7]: 5%  gu
                total_pairs = design['gu'] + design['gc'] + design['ua']
		gcPercentage = float(design['gc'])/total_pairs
		uaPercentage = float(design['ua'])/total_pairs
		guPercentage = float(design['gu'])/total_pairs

		if(gcPercentage>params[3] and gcPercentage<params[4]):
                        score += params[8]
                elif(uaPercentage>params[5] and uaPercentage<params[6]):
                        score += params[8]
                elif(guPercentage<params[7]):
                        score += params[8]
                #print(score)
		return score

