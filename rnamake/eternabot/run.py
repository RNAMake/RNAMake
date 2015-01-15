import sys
import os
import re
import run_utils
import eterna_utils

path="strategies/"
publish =""
op = ""

strategyname = sys.argv[1]



op = ""
if(len(sys.argv) >=3):
	op = sys.argv[2]

print "QUERYING server"
designs = eterna_utils.get_synthesized_designs_from_eterna_server(False, "http://eterna.cmu.edu/eterna_get_synthesized.php?all=1")
print "DONE parsing"

csv_file = open("strategies.csv","w+")
	
if(strategyname == "ALL"):
	published_count = 0
	comprehensive_count = 0
	dirList=os.listdir(path)
	for fname in dirList:
		if(re.search('\.py$',fname)):
			strategy = eterna_utils.load_strategy_from_file('strategies/' + fname)
			if(strategy.publishable_):
				published_count +=1
			if(strategy.author_ == "NUPACK"):
				continue
			run_utils.run_strategy(strategy,designs,op)
	
	print "PUBLISHED : %d" % published_count
else:
	strategy = 	eterna_utils.load_strategy_from_file(strategyname)
	print "Loaded strategy " + strategy.title_
	run_utils.run_strategy(strategy,designs,op)
	
csv_file.close()	