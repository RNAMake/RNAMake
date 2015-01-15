import eterna_utils
import test_utils
import inv_utils
import random
import sys
import re
import os
import math
import numpy
import scipy.linalg
import simplejson
import settings

class Ensemble:
    def __init__(self, method, strategy_names, weights):
        try:
            base_dir = settings.base_dir
            fw = open(base_dir + "/strategies/weights_%s.txt" % method, "r")
            self.my_weights_ = simplejson.loads(fw.read())
            fw.close()

            fs = open(base_dir + "/strategies/strategies_%s.txt" % method, "r")
            strategy_stored_names = simplejson.loads(fs.read())
            fs.close()

            fn = open(base_dir + "/strategies/normalized/mean_%s.txt" % method, "r")
            self.scores_mean_ = simplejson.loads(fn.read())
            fn.close()

            fst = open(base_dir + "/strategies/normalized/stdev_%s.txt" % method, "r")
            self.scores_stdev_ = simplejson.loads(fst.read())
            fst.close()

            if strategy_names != strategy_stored_names:
                raise IOError

            strategies = []

            for i in range(0, len(strategy_stored_names)):
                strategy = eterna_utils.load_strategy_from_file(base_dir + '/strategies/' + strategy_stored_names[i] + ".py")
                strategy.load_opt_params()
                strategy.get_normalization(strategy.default_params_)
                strategies.append(strategy)

            self.strategies_ = strategies

        except IOError:
            exit()
            designs = eterna_utils.get_synthesized_designs_from_eterna_server(True, "http://eterna.cmu.edu/eterna_get_synthesized.php?all=1")
            scores_mean = designs[0]['normalize_mean']
            scores_stdev = designs[0]['normalize_stdev']

            self.designs_ = designs
            self.scores_mean_ = scores_mean
            self.scores_stdev_ = scores_stdev

            strategies = []

            for ii in range(0,len(strategy_names)):
                strategy = eterna_utils.load_strategy_from_file('strategies/' + strategy_names[ii] + ".py")
                strategy.load_opt_params()
                strategy.set_normalization(designs, strategy.default_params_)
                strategies.append(strategy)

            my_weights = []

            # Sparse
            if(weights):
                my_weights_walker = 0
                my_feature_array = []
                my_scores_array = []

                for ii in range(0,len(designs)):

                    my_feature_vector = []
                    for ss in range(0,len(strategies)):
                        if(weights[ss] != 0):
                            f_score =  strategies[ss].normalized_score(designs[ii],strategies[ss].default_params_)['normalized']
                            my_feature_vector.append(f_score)

                    my_feature_vector.append(1.0)
                    my_feature_array.append(my_feature_vector)
                    my_scores_array.append([designs[ii]['score']])

                my_feature_matrix = numpy.matrix(my_feature_array)
                my_scores_matrix = numpy.matrix(my_scores_array)
                my_weights_res = scipy.linalg.lstsq(my_feature_matrix, my_scores_matrix)

                for ss in range(0,len(strategies)+1):
                    if(weights[ss] != 0):
                        my_weights.append(my_weights_res[0][my_weights_walker][0])
                        my_weights_walker += 1
                    else:
                        my_weights.append(0)

                num_of_non_zeros = 0
                for ww in range(0,len(weights)):
                    if(weights[ww] != 0):
                        num_of_non_zeros += 1

                if(my_weights_walker != num_of_non_zeros):
                    print "Error in sparse ensemble %d %d %d" % (my_weights_walker, num_of_non_zeros, len(strategies)+1)
                    sys.exit(0)
            #L2
            else:
                my_feature_array = []
                my_feature_array_t = []
                my_scores_array = []

                for ii in range(0,len(designs)):

                    my_feature_vector = []
                    for ss in range(0,len(strategies)):
                        f_score =  strategies[ss].normalized_score(designs[ii],strategies[ss].default_params_)['normalized']
                        my_feature_vector.append(f_score)

                    my_feature_vector.append(1.0)
                    my_feature_array.append(my_feature_vector)
                    my_scores_array.append([designs[ii]['score']])

                for ss in range(0,len(strategies)+1):
                    my_feature_vector_t = []
                    for ii in range(0,len(designs)):
                        my_feature_vector_t.append(my_feature_array[ii][ss])
                    my_feature_array_t.append(my_feature_vector_t)

                lamb = 85.95482
                lamb_array = []
                for ii in range(0,len(strategies)+1):
                    lamb_vector = []
                    for jj in range(0,len(strategies)+1):
                        lamb_vector.append(0)
                    lamb_vector[ii] = lamb
                    lamb_array.append(lamb_vector)

                my_feature_matrix = numpy.matrix(my_feature_array)
                my_feature_matrix_t = numpy.matrix(my_feature_array_t)
                my_scores_matrix = numpy.matrix(my_scores_array)
                lamb_diagonal = numpy.matrix(lamb_array)

                temp = (my_feature_matrix_t * my_feature_matrix + lamb_diagonal)
                tempi = scipy.linalg.pinv(temp)
                tempi = tempi * my_feature_matrix_t * my_scores_matrix

                for ii in range(0,len(strategies)+1):
                    my_weights.append(float(tempi[ii][0]))

            fw = open("strategies/weights_%s.txt" % method, "w")
            fs = open("strategies/strategies_%s.txt" % method, "w")
            fn = open("strategies/normalized/mean_%s.txt" % method, "w")
            fst = open("strategies/normalized/stdev_%s.txt" % method, "w")

            fw.write(simplejson.dumps(my_weights))
            fs.write(simplejson.dumps(strategy_names))
            fn.write(simplejson.dumps(scores_mean))
            fst.write(simplejson.dumps(scores_stdev))

            fw.close()
            fs.close()
            fn.close()
            fst.close()

            self.my_weights_ = my_weights
            self.strategies_ = strategies

    def test_weights(self,weights):
        my_weights = self.my_weights_
        strategies = self.strategies_
        weights_error = 0
        weights_max_error = 0
        for ss in range(0,len(strategies)+1):
            weights_error += math.fabs(my_weights[ss] - weights[ss])
            if(weights_max_error < math.fabs(my_weights[ss] - weights[ss])):
                weights_max_error = math.fabs(my_weights[ss] - weights[ss])

            if(ss <len(strategies)):
                print "%s, %f" % (strategies[ss].title_, weights[ss])

        print "MY WEIGHTS AVG ERROR %f" % (weights_error / float(len(strategies)+1))
        print "MY WEIGHTS MAX ERROR %f" % (weights_max_error)
        print "MY WEIGHTS %s" % my_weights

    def score(self, design):

        scores_mean = self.scores_mean_
        scores_stdev = self.scores_stdev_
        strategies = self.strategies_
        my_weights = self.my_weights_

        scoresum = 0
        scoremap = {}
        for ss in range(0,len(strategies)):
            if(my_weights[ss] != 0):
                scores = strategies[ss].normalized_score(design,strategies[ss].default_params_)
                score = scores['normalized']
                scoresum += score * my_weights[ss]
                #scoremap[strategies[ss].title_] = score * my_weights[ss];
                scoremap[strategies[ss].title_] =  scores['unnormalized']

        scoresum += my_weights[len(strategies)]
        scoremap['finalscore'] = scoresum * scores_stdev + scores_mean
        scoremap['finalscore_normalized'] = scoresum
        return scoremap

    def printtest(self,str):
        print str

    def test_scores(self, scores):

        designs = self.designs_
        scores_mean = self.scores_mean_
        scores_stdev = self.scores_stdev_

        errorsum = 0
        error_max = 0

        score_over_90 = 0
        score_over_85 = 0
        score_over_90_wrong = 0
        score_over_85_wrong = 0
        score_80s = 0
        score_80s_winner = 0
        for ii in range(0,len(designs)):
            amy_score  = self.score(designs[ii])['finalscore']
            errorsum += math.fabs(scores[ii] - amy_score)

            if(amy_score < 90 and amy_score >= 80):
                score_80s += 1
                if(designs[ii]['score'] * scores_stdev + scores_mean >= 94):
                    score_80s_winner += 1

            if(amy_score >= 90):
                actual_score = designs[ii]['score'] * scores_stdev + scores_mean
                score_over_90 += 1
                if(actual_score < 90):
                    score_over_90_wrong += 1

            if(amy_score >= 85):
                actual_score = designs[ii]['score'] * scores_stdev + scores_mean
                score_over_85 += 1
                if(actual_score < 85):
                    score_over_85_wrong += 1

            if(error_max < math.fabs(scores[ii] - amy_score)):
                error_max = math.fabs(scores[ii] - amy_score)

            if(math.fabs(scores[ii] - amy_score) > 2):
                print "%d %f %f %s" % (ii, scores[ii], amy_score, designs[ii]['soltitle'])

        print "AVG SCORE ERROR %f %f %f" % (errorsum / len(scores), errorsum, len(scores))
        print "MAX SCORE ERROR %f" % (error_max)
        print "OVER 90 ERROR  %d / %d" % (score_over_90_wrong, score_over_90)
        print "OVER 85 ERROR  %d / %d" % (score_over_85_wrong, score_over_85)

        print "80s %d / %d" % (score_80s, score_80s_winner)

