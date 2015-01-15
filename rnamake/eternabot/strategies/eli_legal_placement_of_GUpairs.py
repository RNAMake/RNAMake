from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template
import re

class Strategy(strategy_template.Strategy):
        def __init__(self):
                strategy_template.Strategy.__init__(self)
		self.title_ = "Legal placement of GU-pairs"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_strategy_market_legal_placements_of_gu_pairs"
		self.default_params_ = [2,1,0.5,0.5,-2,-4,-2,-1,2,1,1,-1,-2,-0.5]
		self.code_length_ = 247
		self.publishable_ = True


        '''
        def isExist(self, array, element):
                for ii in range(0,len(array)):
                                if(array[ii]==element):
                                        return True
                return False
	def getNum(self, source):
                result="";
                for ii in range(0, len(source)):
                        result= result+source[ii];
                return int(result);
        def is_same_neck(self, elements, pairIndexA, pairIndexB):
                #print("SameNeck",pairIndexA,pairIndexB);
                for ii in range(0,len(elements)):
                        elem = elements[ii]
                        if(elem.type_ == RNAELEMENT_STACK):
                                for ii in range(0,len(elem.indices_)):
                                        if(elem.indices_[ii]==pairIndexA):# find the pairIndexA's position
                                                #print("A's stack:",elem.indices_,"Search:",pairIndexB, "Length:", elem.get_stack_length());
                                                for jj in range(0,len(elem.indices_)):
                                                        if(elem.indices_[jj]==pairIndexB):
                                                                return True;

                return False;
        def get_same_stack(self, elements, pairIndex):
                for ii in range(0,len(elements)):
                        elem = elements[ii]
                        if(elem.type_ == RNAELEMENT_STACK):
                                for ii in range(0,len(elem.indices_)):
                                        if(elem.indices_[ii]==pairIndex):
                                                return elem;

                return 0;

        '''

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
		score = 80
                '''

                indexOfNeck=-1;
                first_nuc_with_pair=-1;
                last_nuc_with_pair=-1;

                for i in range(len(sequence)):
                        if(i!=-1):
                                first_nuc_with_pair=i;
                                last_nuc_with_pair=pairmap[i];

                neckArea= 0;
                for ii in range(0,len(elements)):
                        elem = elements[ii]
                        if(elem.type_ == RNAELEMENT_STACK):
                                for jj in range(0,elem.get_stack_length()):
                                        if(elem.indices_[jj]==first_nuc_with_pair):
                                                break;
                                neckArea=elem;
                                break;
                '''

                noGU=True;
                for ii in range(0,len(elements)):
                        elem = elements[ii]
                        if(elem.type_ == RNAELEMENT_STACK):
                                if(len(self.get_original_pairs(elem, sequence))!=0):
                                        noGU=False;
                                        break;
                if(noGU):
                        score += params[0];

                for ii in range(0,len(elements)):
                        elem = elements[ii]
                        if(elem.type_ == RNAELEMENT_STACK):
                                numBetween=self.two_GC(elem, sequence);
                                score += numBetween[0]*params[1];
                                score += numBetween[1]*params[2];
                                numBeside = self.besideGC(elem, sequence);
                                if(numBeside!=0):
                                        score += numBeside*params[3];
                                else:
                                        score += params[4];

                besidingGU=0;
                for ii in range(0,len(elements)):
                        elem = elements[ii]
                        if(elem.type_ == RNAELEMENT_STACK):
                                besidingGU += self.get_besidingGU_pairs(elem, sequence);
                if(besidingGU==1):
                        score += params[5];
                elif(besidingGU>1):
                        score += (besidingGU-1)*params[6];


                for ii in range(0,len(elements)):
                        elem = elements[ii];
                        if(elem.type_ == RNAELEMENT_STACK):
                                pairs = self.get_original_pairs(elem, sequence);
                                if(len(pairs)==0):
                                        continue;
                                if(len(pairs[0])==0 or len(pairs[len(pairs)-1])==0):
                                        score += params[7];





                indexOfNeck=-1;
                first_nuc_with_pair=-1;
                last_nuc_with_pair=-1;

                for i in range(len(sequence)):
                        if(i!=-1):
                                first_nuc_with_pair=i;
                                last_nuc_with_pair=pairmap[i];

                neckArea= 0;
                for ii in range(0,len(elements)):
                        elem = elements[ii]
                        if(elem.type_ == RNAELEMENT_STACK):
                                for jj in range(0,elem.get_stack_length()):
                                        if(elem.indices_[jj]==first_nuc_with_pair):
                                                break;
                                neckArea=elem;
                                break;




                numBetween=self.two_GC(neckArea, sequence);
                score += numBetween[0]*params[8];
                score += numBetween[1]*params[9];
                numBeside = self.besideGC(neckArea, sequence);
                if(numBeside!=0):
                        score += numBeside*params[10];
                else:
                        score += params[11];


                pairs=[];
                tmp=0;
                isSameTuning=True;
                for i in range(0,neckArea.get_stack_length()):
                        pair = neckArea.get_pair_from_stack(i,sequence);
                        pairs.append(pair);
                for i in range(0,len(pairs)-1):
                        if((pairs[i]=="GU" or pairs[i]=="UG") and\
                           (pairs[i+1]=="GU" or pairs[i+1]=="UG")):
                                tmp += 1;
                                if(pairs[i]!=pairs[i+1]):
                                        isSameTurning=False;
                if(tmp >1 or isSameTuning):
                        score += params[12];

                if(pairs[0]=="GU" or pairs[0]=="UG" or pairs[len(pairs)-1]=="GU" or pairs[len(pairs)-1]=="UG"):
                        score += params[13];

                #print(neckArea.get_stack_length());
                #self.two_GC(elements[0], sequence);
                #print(self.besideGC(elements[0], sequence));
                #print(self.get_besidingGU_pairs(elements[0], sequence));
                return score



        def get_num_twistGC(self, groupA,groupB):
                result=0;
                for pairA in groupA:
                        for pairB in groupB:
                                if(pairA!=pairB):
                                        result = result+1;
                return result;
        def get_num_sameGC(self, groupA,groupB):
                result=0;
                for pairA in groupA:
                        for pairB in groupB:
                                if(pairA==pairB):
                                        result = result+1;
                return result;


        def get_besidingGU_pairs(self, elem, sequence):
                pairs=[];
                result = 0;
                for i in range(0,elem.get_stack_length()):
                        pair = elem.get_pair_from_stack(i,sequence);
                        pairs.append(pair);
                for i in range(0,len(pairs)-1):
                        if((pairs[i]=="GU" or pairs[i]=="UG") and\
                           (pairs[i+1]=="GU" or pairs[i+1]=="UG")):
                                result += 1;
                return result

        def get_original_pairs(self, elem, sequence):
                pairs=[];
                pairs_sep=[];
                for i in range(0,elem.get_stack_length()):
                        pair = elem.get_pair_from_stack(i,sequence);
                        pairs.append(pair);
                tmp=0;
                p=[];
                for i,pair in enumerate(pairs):
                        if(pair=="GU" or pair=="UG"):
                                p = [];
                                for j in range(tmp,i):
                                        p.append(pairs[j]);

                                if(len(p)>0):
                                        pairs_sep.append(p);
                                elif(tmp==0):
                                        pairs_sep.append(p);

                                tmp=i+1;
                p=[];
                if(tmp!=0):
                        for j in range(tmp,len(pairs)):
                                p.append(pairs[j]);
                        pairs_sep.append(p);

                #print("\n");
                #print("================");
                #print(pairs);

                return pairs_sep;

        def get_filted_pairs(self, elem, sequence):

                pairs_sep=self.get_original_pairs(elem, sequence);
                #print(pairs_sep);


                pairs_sep_filted=[];
                for ii in range(len(pairs_sep)):
                        tmp=[];
                        for jj in range(len(pairs_sep[ii])):
                                if(pairs_sep[ii][jj]=="GC" or pairs_sep[ii][jj]=="CG"):
                                        tmp.append(pairs_sep[ii][jj]);
                        pairs_sep_filted.append(tmp);
                #print(pairs_sep_filted);
                return pairs_sep_filted;

        def besideGC(self, elem, sequence):#return number of besiding one GCpair and no GCpair
                result = 0;
                pairs= self.get_original_pairs(elem, sequence);

                if(len(pairs)<2):
                        return result;

                tmp= len(pairs[0]);
                if(tmp!=0 and (pairs[0][tmp-1]=="GC" or pairs[0][tmp-1]=="CG")):
                        #print("Y",pairs[0][tmp-1]);
                        result += 1;

                if(len(pairs)>0):
                        tmp= len(pairs)-1;
                        if(len(pairs[tmp])>0  and (pairs[tmp][0]=="GC" or pairs[tmp][0]=="CG")):
                                #print("N",pairs[tmp][0]);
                                result += 1;

                if(len(pairs)>2):
                        for ii in range(1, len(pairs)-1):
                                tmp= len(pairs[ii]);

                                if(tmp>0 and (pairs[ii][tmp-1]=="GC" or pairs[ii][tmp-1]=="CG")):
                                        #print("III", pairs[ii][tmp-1]);
                                        result += 1;

                                if(len(pairs[ii])>0 and(pairs[ii][0]=="GC" or pairs[ii][0]=="CG")):
                                        #print("---", pairs[ii][0]);
                                        result += 1;
                #print("??", result)
                return result;




        def two_GC(self, elem, sequence):#return the number of twisted and same turning pairs

                numOfTwist=0;
                numOfSame=0;
                result=[];
                pairs_sep_filted = self.get_filted_pairs(elem, sequence);
                for ii in range(0,len(pairs_sep_filted)-1):
                        pair = pairs_sep_filted[ii];
                        for jj in range(ii+1, len(pairs_sep_filted)):
                                compare_pair=pairs_sep_filted[jj];
                                #print(pair, compare_pair);
                                numOfTwist += self.get_num_twistGC(pair,compare_pair);
                                numOfSame += self.get_num_sameGC(pair,compare_pair);
                result.append(numOfTwist);
                result.append(numOfSame);

                #print(result);
                return result;






