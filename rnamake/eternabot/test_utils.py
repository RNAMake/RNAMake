import re
import random
import math
import eterna_utils
import os
from scipy.optimize import fmin, fmin_powell
#from scikits.learn import svm


def optimize_strategy_with_squared_error(strategy, training_set): 
	if(len(strategy.default_params_) == 0):
		return strategy.default_params_
	
	def score_func(params): 
		return calculate_average_squared_error(params,training_set,strategy)
	
	return fmin(score_func, strategy.default_params_)
	
	
def optimize_strategy_with_kt_corr(strategy, training_set): 
	if(len(strategy.default_params_) == 0):
		return strategy.default_params_
	
	def score_func(params): 
		return -calculate_kt_corr(params,training_set,strategy)
	
	return fmin(score_func, strategy.default_params_)
	
	

def cross_validate(n,designs,strategy):
	num_designs = len(designs)
	test_num = int(num_designs / float(n))
	
	score_sum = 0.0
	
	for ii in range(0,n):
		
		test_set = []
		training_set = []
		
		for jj in range(0,num_designs):
			if(jj >= test_num * ii and  jj < test_num * (ii+1)):
				test_set.append(designs[jj])
			else:
				training_set.append(designs[jj]);

		opt_params = optimize_strategy_with_kt_corr(strategy,training_set)
		this_score = calculate_kt_corr(opt_params,test_set, strategy)
		score_sum += float(this_score)
	
	return score_sum / float(n)

def cross_validate_score_tuples(n,designs,strategy):
	num_designs = len(designs)
	num_tests = int(num_designs / float(n))
	
	test_set = []
	training_set = []
	
	score_tuples = []
	
	for ii in range(0,num_designs):
		if(ii < num_designs - num_tests):
			training_set.append(designs[ii])
		else:
			test_set.append(designs[ii])
	
	opt_params = optimize_strategy_with_kt_corr(strategy,training_set)
	
	for ii in range(0,num_designs):
		if(ii < num_designs - num_tests):
			continue
		else:
			prediction = (strategy.score(designs[ii], opt_params))
			score_tuple = {}
			score_tuple['score'] = designs[ii]['score']
			score_tuple['predicted_score'] = prediction
			score_tuples.append(score_tuple)
	
	return score_tuples


def leave_one_out_cv_classifier(params, designs,strategies):
	puznids = get_puz_nids(designs)
	outs = []
	sum = 0
	
	if(len(strategies) > 1):
		title = "Ensemble %d" % len(strategies)
	else:
		title= strategies[0].title_
	
	for ii in range(0,len(puznids)):
		set = get_leave_one_out_set(designs,puznids[ii])
		cl = get_classifier(set['train'],strategies,params)
		err = calculate_classifier_error(params,set['test'],strategies,cl)
		sum += err
		print "%s PUZ %d :: %f" % (title, puznids[ii], err)
		outs.append(err)
	
	
	res = {}
	res['avg'] = sum /len(puznids)
	res['list'] = outs
	res['title'] = title
	
	print " "
	
	return res	

def get_classifier(designs, strategies, params):
	score_vectors = []
	res_vector = []
	
	slm = success_label_multiplier(designs)
	
	for ii in range(0,len(designs)):
		score_vector = []
		for jj in range(0,len(strategies)):
			score_vector.append(strategies[jj].score(designs[ii],params[jj]))
		
		score_vectors.append(score_vector)
		if(designs[ii]['score'] >= 90):
			res_vector.append(1)
		else:
			res_vector.append(-1)
			
	#clf = svm.SVC()
	clf.fit(score_vectors, res_vector, class_weight={1:slm,-1:1}, kernel="linear")

	return clf

def calculate_classifier_error(params, input_designs,strategies, classifier):
	i_range = range(0,len(input_designs)) 
	j_range = range(0,len(strategies))

	error_count = 0
	total_count = 0

	slm = success_label_multiplier(input_designs)

	for ii in i_range:
		score_vector = []
		for jj in j_range:
			score_vector.append(strategies[jj].score(input_designs[ii],params[jj]))
		
		#print str(score_vector)
		
		res = classifier.predict([score_vector])

		mul = 1
		if(input_designs[ii]['score'] >= 90):
			mul = slm
		
		total_count += 1 * mul
		
		if(res < 0 and input_designs[ii]['score'] >= 90):
			error_count += 1 * mul
		if(res > 0 and input_designs[ii]['score'] < 90):
			error_count += 1 * mul
			
	
	return error_count / float(total_count)
	

def calculate_average_squared_error(params, input_designs,strategy): 
	i_range = range(0,len(input_designs)) 
	error_sum = 0
	count = 0

	for ii in i_range:
		prediction = (strategy.score(input_designs[ii], params))
		synthesis = input_designs[ii]['score']
		input_designs[ii]['predicted_score'] = prediction
	
		if(prediction <= eterna_utils.UNSCORABLE):
			continue
		
		count += 1
		error_sum += (prediction - synthesis) * (prediction - synthesis)
		
	
	if(count > 0):
		error_sum = float(error_sum)/count
		
	return error_sum

def calculate_average_plot_error(params, input_designs, strategy):
	i_range = range(0,len(input_designs)) 
	error_sum = 0
	count = 0

	for ii in i_range:
		prediction = (strategy.score(input_designs[ii], params))
		synthesis = input_designs[ii]['score']
		input_designs[ii]['predicted_score'] = prediction
	
		if(prediction <= eterna_utils.UNSCORABLE):
			continue
		
		count += 1
		error_sum += (math.floor(prediction/10.0) - math.floor(synthesis/10.0)) * (math.floor(prediction/10.0) - math.floor(synthesis/10.0)) 
	
	if(count > 0):
		error_sum = float(error_sum)/count
		
	return error_sum
	
def calculate_kt_corr(params, input_designs, strategy):

	i_range = range(0,len(input_designs))
	wrong_count = 0
	correct_count = 0
	total_count =0;

	for ii in i_range:
		input_designs[ii]['predicted_score'] = (strategy.score(input_designs[ii], params))

	for ii in i_range:
		j_range = range(ii+1, len(input_designs))
		for jj in j_range:
			if (input_designs[ii]['score'] < input_designs[jj]['score']):
				actual_comp = 1
			elif (input_designs[ii]['score'] > input_designs[jj]['score']):
				actual_comp = -1
			else:
				actual_comp = 0
				
			if(math.fabs(input_designs[ii]['score'] - input_designs[jj]['score']) <= 4.0):
				actual_comp = 0

			if (input_designs[ii]['predicted_score'] < input_designs[jj]['predicted_score']):
				pred_comp = 1
			elif (input_designs[ii]['predicted_score'] > input_designs[jj]['predicted_score']):
				pred_comp = -1
			else:
				pred_comp = 0
				
			if(input_designs[ii]['predicted_score'] <= eterna_utils.UNSCORABLE or input_designs[jj]['predicted_score'] <= eterna_utils.UNSCORABLE):
				pred_comp = 0
			
			if(pred_comp * actual_comp != 0):
				if(pred_comp != actual_comp):
					wrong_count +=1
				else:
					correct_count +=1
			
			total_count += 1

	kt_corr = float(correct_count - wrong_count) / total_count
	return kt_corr




def get_score_plot_code(designs,strategy):
	ret =""
	
	ret += ('<script language="javascript" type="text/javascript" src="http://eterna.cmu.edu/jss/flot/jquery.js"></script>\n')
	ret += ('<script language="javascript" type="text/javascript" src="http://eterna.cmu.edu/jss/flot/jquery.flot.js"></script>\n\n')

	plotid = 'plot' + str(int(random.random() * 1000))

	ret += ("<div id='" + plotid + "' style='float:left; width:400px; height:400px; margin-top:20px'></div>")
	ret += ("<script type='text/javascript'>\n")
	ret += ("$(function () {\n\n")
	ret += ("var d1 = [ ")
	for ii in range(0, len(designs)):
		if ii > 0:
			ret += (", ")
	
		if(designs[ii]['predicted_score'] > eterna_utils.UNSCORABLE and designs[ii]['score'] > 20): 
			ret += ("[" + str(designs[ii]['score']) + "," + str(designs[ii]['predicted_score']) + "]")

		
	ret += ("];\n")
	ret += ("$.plot($('#" + plotid + "'), [")
	ret += ("{ label: '" + strategy.title_ + "', data:d1, points : {show:true, radius:5}, color:'#00F' }")
	ret += ("],")
	ret += ("{ xaxis : {max: 100, min: 0}, yaxis : {max: 100, min: 0},grid: { backgroundColor: { colors: ['#fff', '#eee'] } } }")
	ret += (");\n\n")
	ret += ("});")

	ret += ("</script>");

	return ret
	
def success_label_multiplier(designs):
	total = 0
	success = 0
	for ii in range(0,len(designs)):
		if(designs[ii]['score'] >= 90):
			success += 1
		total += 1

	return float(total) / success
	
	
def get_puz_nids(designs):
	nids = []
	for ii in range(0,len(designs)):
		if(designs[ii]['puznid'] not in nids):
			nids.append(designs[ii]['puznid'])
	
	return nids
	
def get_leave_one_out_set(designs,puznid):
	train = []
	test = []
	
	for ii in range(0,len(designs)):
		if(designs[ii]['puznid'] == puznid):
			test.append(designs[ii])
		else:
			train.append(designs[ii])
	res = {}
	res['test'] = test
	res['train'] = train
	return res




