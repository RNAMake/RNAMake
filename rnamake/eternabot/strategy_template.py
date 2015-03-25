import eterna_utils
import test_utils
import sys
import math
import json
import simplejson
import settings

class Strategy(json.JSONEncoder):
    def __init__(self):
        self.title_ = ""
        self.author_ = ""
        self.url_ = ""
        self.default_params_ = []
        self.code_length_ = 0
        self.publishable_ = False
        self.is_part_of_ensemble_ = False
        self.optimized_ = False
        self.denormalized_ = False
        self.comprehensive_ = True
        self.martin_weight_ = 0
        self.satisfying_point_ = None

        self.mean_ = None
        self.stdev_ = None

    def load_opt_params(self):

        params_name = settings.base_dir + "/strategies/params/" + self.author_ + "_" + self.title_ + ".params"
        try :
            params_file = open(params_name,"r")
            params_content = params_file.read()
            params_file.close()
            params = params_content.split()

            if(len(params) != len(self.default_params_)):
                print "Failed to load params - # of parameters in " + params_name + " does not match that of default_params"
                return

            for ii in range(0,len(params)):
                self.default_params_[ii] = float(params[ii])

            self.optimized_ = True
            return True

        except IOError :
            print "Failed to load params - params file does not exist"
            return False

    def save_params(self, opt_params):
        params_name = "./strategies/params/" + self.author_ + "_" + self.title_ + ".params"
        try :

            if(len(opt_params) != len(self.default_params_)):
                print "Failed to save params - # of parameters in opt_params does not match that of default_params"
                return

            params_file = open(params_name,"w+")



            for ii in range(0,len(opt_params)):
                if(ii>0):
                    params_file.write(" ")
                params_file.write(str(opt_params[ii]))

            params_file.close()

            return

        except IOError :
            print "Failed to save params - could not open file to write somehow.."

    def save_cv_kt_score(self, n, ktscore):
        score_name = "./strategies/scores/" + self.author_ + "_" + self.title_ + ".cvkt" + str(n)

        try :

            score_file = open(score_name,"w+")
            score_file.write(str(ktscore))
            score_file.close()

            self.cvkt = ktscore

        except IOError :
            print "Failed to save kt score - could not open file to write somehow.."

    def load_cv_kt_score(self, n):
        score_name = "./strategies/scores/" + self.author_ + "_" + self.title_ + ".cvkt" + str(n)

        try :
            score_file = open(score_name,"r")
            score_content = score_file.read()
            self.cvkt = float(score_content)
            return True
        except (IOError, ValueError) :
            print "Failed to load kt score - score file does not exist"
            return False

    def get_normalization(self, params):

        strategy_name = self.author_ + "_" + self.title_

        fm = open(settings.base_dir + "/strategies/normalized/" + strategy_name + "_mean.txt", "r");
        fs = open(settings.base_dir + "/strategies/normalized/" + strategy_name + "_stdev.txt", "r");

        self.mean_ = simplejson.loads(fm.read())
        self.stdev_ = simplejson.loads(fs.read())

        fm.close()
        fs.close()

    def set_normalization(self, designs, params):
        scores = []
        scoresum = 0
        n = 0

        for ii in range(0,len(designs)):
            s = (self.score(designs[ii], params))
            if(s <= eterna_utils.UNSCORABLE):
                s = 0
            scoresum += s
            n += 1
            scores.append(s)

        mean = scoresum
        if(n > 0):
            mean = mean / float(n)

        stdev = 0
        for ii in range(0,len(designs)):
            stdev += (scores[ii] - mean) * (scores[ii] - mean)

        if(n > 0):
            stdev = stdev / float(n)

        stdev = math.sqrt(stdev)

        strategy_name = self.author_ + "_" + self.title_

        fm = open("strategies/normalized/" + strategy_name + "_mean.txt", "w")
        fs = open("strategies/normalized/" + strategy_name + "_stdev.txt", "w")

        fm.write(simplejson.dumps(mean))
        fs.write(simplejson.dumps(stdev))

        fm.close()
        fs.close()

        self.mean_ = mean
        self.stdev_ = stdev

        #print "MS %s %f %f" % (self.title_, self.mean_, self.stdev_)

    def normalized_score(self, design,params):
        if(self.mean_ == None or self.stdev_ == None):
            print "Strategy not normalized"
            sys.exit(0)

        s = (self.score(design, params))
        orig_s = s
        if(s <= eterna_utils.UNSCORABLE):
            s = 0

        #print self.title_, self.mean_, self.stdev_
        return {"normalized": (s - self.mean_) / (self.stdev_ + 1e-15), "unnormalized": (orig_s)}


    def score(self, design, params):
        return 0

    def patch(self, design, params):
        return None

    def get_title(self, obj):
        return self.title_
