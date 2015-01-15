strategies = []
strategies.append(eterna_utils.load_strategy_from_file("strategies\\xmbrst_clear_plot_stack_caps_and_safe_gc.py"))
#strategies.append(eterna_utils.load_strategy_from_file("strategies\\aldo_repetition.py"))
strategies.append(eterna_utils.load_strategy_from_file("strategies\\ding_quad_energy.py"))
strategies.append(eterna_utils.load_strategy_from_file("strategies\\deivad_deivad_strategy.py"))

"""
path="strategies/"
dirList=os.listdir(path)
strategies = []

for fname in dirList:
	if(re.search('\.py$',fname) and re.search('ensemble',fname) == None and re.search('nupack',fname) == None):
		strategy = eterna_utils.load_strategy_from_file('strategies\\' + fname)
		if(strategy.denormalized_):
			strategies.append(strategy)
"""


params = []
for ii in range(0,len(strategies)):
	params.append(strategies[ii].default_params_)

designs = eterna_utils.get_synthesized_designs_from_eterna_server()


	
"""
results = []

for ii in range(0,len(strategies)):
	results.append(leave_one_out_cv_classifier([params[ii]],designs,[strategies[ii]]))

results = sorted(results, key=lambda st: st['avg'], reverse=True)

for ii in range(0,len(results)):
	print "%s %f %s" % (results[ii]['title'], results[ii]['avg'], str(results[ii]['list']))
"""

print str(leave_one_out_cv_classifier(params, designs, strategies))