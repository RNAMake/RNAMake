import sys
import re
import urllib
import random
import eterna_utils
import test_utils
import httplib2
import json
from scipy.optimize import fmin


def run_strategy(strategy,designs,op):

	if(op != "OPTIMIZE"):
		strategy.load_opt_params()
	else:
		print "Redoing optimization.."
		
	print "Testing " + strategy.author_ + "'s " + strategy.title_
	print "OPTIMIZED " + str(strategy.optimized_)
	
	
	file_prefix = ""
	f_output = open("./strategies/outputs/" + file_prefix + strategy.author_ + "'s " + strategy.title_ + ".html","wb")
		
	
	if(op == "CROSS"):
		CROSS_N = 3		
		return test_utils.cross_validate(CROSS_N,designs,strategy)

	
	if(len(strategy.default_params_) > 0 and (strategy.optimized_ == False or op == "OPTIMIZE")):
		#opt_params = test_utils.optimize_strategy_with_kt_corr(strategy,designs)
		opt_params = test_utils.optimize_strategy_with_squared_error(strategy,designs)
		final_kt_corr = test_utils.calculate_kt_corr(opt_params, designs, strategy)
		default_kt_corr = test_utils.calculate_kt_corr(strategy.default_params_, designs, strategy)
	else:
		opt_params = strategy.default_params_
		final_kt_corr = test_utils.calculate_kt_corr(strategy.default_params_, designs, strategy)
		default_kt_corr = final_kt_corr
	
	print "KT correlation calculated"
	
	strategy.save_params(opt_params)
	
	 
	
	print("FINAL " + str(final_kt_corr))
	
	
	
	## Output results
	f_output.write("<strong>Title : </strong>" + strategy.title_ + "<br/>")
	f_output.write("<strong>Author : </strong>" + strategy.author_ + "<br/>")
	f_output.write("<strong>Final KT corr : </strong>" + str(final_kt_corr) + "<br/>")
	f_output.write("<strong>Final params : </strong>" + str(opt_params) + "<br/>")
	f_output.write("<strong>Initial KT corr : </strong>" + str(default_kt_corr) + "<br/>")
	f_output.write("<strong>Initial params : </strong>" + str(strategy.default_params_) + "<br/><br/>")
	
	f_output.write("\n\n\n")
	
	body_data = ""
	
	body_data += ("<table style='border-style:none; padding-top:0px; margin:0px; width:100%'>")
	body_data += ("<tr style='padding-bottom:11px'>")
	body_data += ("<th>Lab name</th>")
	body_data += ("<th>Design name</th>")
	body_data += ("<th>Predicted score</th>")
	body_data += ("<th>Synthesis score</th>")
	body_data += ("<th>Error</th>")
	body_data += ("</tr>")
	
	i_range = range(0,len(designs))
	for ii in i_range: 
		body_data += ("<tr>")
		body_data += ("<td>")
		body_data += ("<a href='http://eterna.cmu.edu/eterna_page.php?page=game&type=BROWSE&value=" + str(designs[ii]['puznid']) +"&filter1=Synthesized&filter1_arg1=y' target='_blank'>")
		body_data += (designs[ii]['puztitle'])
		body_data += ("</a>")
		body_data += ("</td>")
		body_data += ("<td>")
		body_data += ("<a href='http://eterna.cmu.edu/eterna_page.php?page=game&type=BROWSE&value=" + str(designs[ii]['puznid']) +"&filter1=Synthesized&filter1_arg1=y&filter2=Id&filter2_arg1=" + str(designs[ii]['solnid']) +"&filter2_arg2=" + str(designs[ii]['solnid']) + "' target='_blank'>")
		body_data += (designs[ii]['soltitle'] + " (id:" + str(designs[ii]['solnid']) + ")")
		body_data += ("</a>")
		body_data += ("</td>")
		body_data += ("<td style='width:50px'>")
		if(designs[ii]['predicted_score'] > eterna_utils.UNSCORABLE):
			body_data += (str(designs[ii]['predicted_score']))
		else:
			body_data += "N/A"
		body_data += ("</td>")
		body_data += ("<td style='width:50px'>")
		body_data += (str(designs[ii]['score']))
		body_data += ("</td>")
		body_data += ("<td style='width:50px; color:#FFAAAA'>")
		if(designs[ii]['predicted_score'] > eterna_utils.UNSCORABLE):
			diff = abs(designs[ii]['predicted_score'] - designs[ii]['score'])
			body_data += (str(diff))
		else:
			body_data += "N/A"
		body_data += ("</td>")
		body_data += ("</tr>")
	
	body_data += ("</table>")
	
	f_output.write(body_data)
	
	
	data = ""
	if(op == "PUBLISH" and strategy.publishable_):
		conn = httplib2.Http(".cache")
		headers = {'Content-type': 'application/x-www-form-urlencoded'}
		body = {'title':strategy.title_, 'author':strategy.author_, 'url':strategy.url_, 'ktcorr':'%.3f' % final_kt_corr, 'body':body_data, 'codelen':strategy.code_length_ , 'type':"strategy",  'workbranch':"main"}
		response, data = conn.request("http://eterna.cmu.edu/eterna_post.php", "POST", headers=headers, body=urllib.urlencode(body))
		print data
	
	
	f_output.write("\n\n\n");
	f_output.write("Dear " + strategy.author_ + "<br/><br/>");
	f_output.write("We are glad to report that your strategy has been implemented and tested.<br/><br/>");
	f_output.write("While implementing your strategy, we have made small changes to the parameters you specified to optimize the performance.<br/><br/>");
	f_output.write("Note that we'll always run a optimization over the parameters you specify, so you won't have to worry about fine tuning all the numbers you use.<br/><br/>");
	f_output.write("Just the idea and rough numbers are enough to run your algorithm!<br/><br/>");
	f_output.write("<strong>Length</strong> : Your strategy was implmented with <strong>" + str(strategy.code_length_) + "</strong> line of code.<br/><br/>");
	f_output.write("<strong>Ordering</strong> : We ran your strategy on all synthesized designs and ordered them based on predicted scores. The correlation of your strategy's ordering with the ordering based on the actual scores was <strong>" + str(final_kt_corr) + "</strong>. (1.0 is the best score, -1.0 is the worst score. A completely random prediction would have 0 correlation)<br/><br/>"); 
	f_output.write("Please note that the numbers specified above will change in future as we'll rerun your algorithm whenever new synthesis data is available.<br/><br/>");
	if(len(data) > 0):
		ret = json.loads(data)
		finalnid = ret['success']
		f_output.write("More detailed result has been posted on the <a href='http://eterna.cmu.edu/eterna_page.php?page=strategy&nid=" + finalnid + "' target='_blank'>strategy market page</a>. Thank you for sharing your idea, and we look forward to other brilliant strategies from you!");
	
	f_output.write("<br/>");
	
	"""
	f_output.write(test_utils.get_score_plot_code(designs,strategy))
	cross_score_tuples = test_utils.cross_validate_score_tuples(3,designs,strategy)
	f_output.write(test_utils.get_score_plot_code(cross_score_tuples,strategy))
	"""
	
	f_output.close();
	
	

