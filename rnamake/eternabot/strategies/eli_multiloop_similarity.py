from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
	def __init__(self):
		strategy_template.Strategy.__init__(self)

		self.title_ = "Multiloop similarity"
		self.author_ = "Eli Fisker"
		self.url_ = "http://getsatisfaction.com/eternagame/topics/_market_strategy_multiloop_similarity"
		self.default_params_ = [0.5,-1]
		self.code_length_ = 14
		self.publishable_ = False
        def isExist(self, array, element):
                for ii in range(0,len(array)):
                                if(array[ii]==element):
                                        return True
                return False
        def getChildInfo(self,elem,out, energy):

                out.append(str(elem.type_)+str(len(elem.indices_)));
                energy.append(str(elem.score_))
                for ii in range(0,len(elem.children_)):
                        self.getChildInfo(elem.children_[ii], out,energy);

	def score(self, design, params):
		elements = design['secstruct_elements']
		score = 100
                '''
		for ii in range(0,len(elements)):

                        elem = elements[ii]

                        if(elem.type_ == RNAELEMENT_LOOP):
                                loop_groups = elem.get_loop_groups()
                                if(len(loop_groups)>=3):
                                        for ii in range(0,len(loop_groups)-1):
                print("#############################");
                elem = elements[0]

                print(elem.type_);
                for ii in range(0,len(elem.children_)):
                        print(elem.type_);
                '''
                #print("#############################");
                root = elements[len(elements)-1]
                infoGrp=[]
                energyGrp=[]
                for ii in range(0, len(root.children_)):
                        #print(root.children_[ii].type_)
                        parent=root.children_[ii].children_[0]
                        loop_groups=parent.get_loop_groups();

                        for jj in range(0,len(parent.children_)):
                                info=[]
                                energy=[]
                                self.getChildInfo(parent.children_[jj],info,energy)
                                infoGrp.append(info)
                                energyGrp.append(energy)
                        if(len(loop_groups)==0):
                                continue;
                        for jj in range(0,len(loop_groups)-1):
                                if(len(loop_groups[jj])!=len(loop_groups[jj+1])):
                                        #print(loop_groups);
                                        return score

                #print(infoGrp)
                identicalGrp=[]
                for ii in range(0, len(infoGrp)):
                        for jj in range(0,len(infoGrp)):
                                if(ii!=jj and infoGrp[ii]==infoGrp[jj]):
                                        tmp=[]
                                        tmp.append(min(ii,jj));
                                        tmp.append(max(ii,jj));
                                        if(self.isExist(identicalGrp,tmp)==False):
                                                identicalGrp.append(tmp);
                #print("IGrp",identicalGrp)
                #print("EGrp",energyGrp)
                penalty=0;
                for ii in range(0, len(identicalGrp)):
                        e1= float(0)
                        for jj in range(0,len(energyGrp[identicalGrp[ii][0]])):
                                e1 = e1 + float(energyGrp[identicalGrp[ii][0]][jj])
                        e2= float(0)
                        for jj in range(0,len(energyGrp[identicalGrp[ii][1]])):
                                e2 = e2 + float(energyGrp[identicalGrp[ii][1]][jj])
                        #print(e1,e2)
                        if(e1!=e2):
                                e=abs(e1-e2)
                                penalty= int(penalty)+int((e-(e%params[0]))/params[0]);
                                #print((e-(e%params[0]))/params[0], int((e-(e%params[0]))/params[0]));

                #print(penalty)
                score = score + penalty* params[1]
                #print(score,",",params[1])
                return score
















