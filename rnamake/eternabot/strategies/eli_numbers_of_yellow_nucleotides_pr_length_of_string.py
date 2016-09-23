from rnamake.eternabot.eterna_utils import *
import rnamake.eternabot.inv_utils as inv_utils
import rnamake.eternabot.strategy_template as strategy_template

class Strategy(strategy_template.Strategy):
    def __init__(self):
        strategy_template.Strategy.__init__(self)

        self.title_ = "Numbers of yellow nucleotides pr length of string"
        self.author_ = "Eli Fisker"
        self.url_ = "http://getsatisfaction.com/eternagame/topics/numbers_of_yellow_nucleotides_pr_length_of_string"
        #ratio of gc pairs in junctions
        self.default_params_ = [200]
        self.upper_length_ = [0,1,2,2,2,3,3,4,5,4]
        self.lower_length_ = [0,0,0,1,1,1,2,3,2,1]
        self.code_length_ = 30
        self.publishable_ = True
        self.denormalized_ = True
        self.comprehensive_ = False

    def score(self, design, params):
        sequence = design['sequence']
        elements = design['secstruct_elements']
        pairmap = design['pairmap']

        penalty= 0
        count = 0
        for ii in range(0,len(elements)):
            if(elements[ii].type_ == RNAELEMENT_STACK and
                elements[ii].get_stack_length()>2 and
                elements[ii].get_stack_length()<10):
                stack_len = elements[ii].get_stack_length()
                count += 1
                indices = elements[ii].indices_;
                yellow_cnt=0
                for jj in range(0, len(indices)):
                    if(sequence[indices[jj]]=="A"):
                        yellow_cnt += 1
                if(self.upper_length_[stack_len]<yellow_cnt):
                    penalty += yellow_cnt - self.upper_length_[stack_len]
                elif(self.lower_length_[stack_len]>yellow_cnt):
                    penalty += self.lower_length_[stack_len] - yellow_cnt
        if(count == 0):
            return 100


        return 100 - params[0] * penalty/len(sequence)

    def patch(self, design, params):


        native = inv_utils.fold(sequence)[0]
        # If too many A, change AUs to GCs, starting from outside
        # If too few A, change GUs then GCs to AU, starting from inside
        return design
