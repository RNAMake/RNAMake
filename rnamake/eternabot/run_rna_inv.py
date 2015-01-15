import csv
import httplib2
import sys
import re
import urllib
import random
import eterna_utils
import inv_utils
import time
#run strategies to inverse fold a structure
MAX_NUM_ROWS = 10
op = "";
if(len(sys.argv) >= 3):
        op = sys.argv[2];


sys.path.append("./rna_inv_strategies")

strategyname = sys.argv[1]
strategyname = strategyname.replace('\\', '/')
m = re.search("[^/]*\.py$", strategyname)
strategyname = strategyname[m.start():m.end() - 3]

import_string = "from " + strategyname + " import *"
exec import_string

strategy = Strategy();

print "\n\n\n============="
print "Testing " + strategy.author_ + "'s " + strategy.title_


file_prefix = ""
f_output = open("./rna_inv_strategies/" + file_prefix + strategy.author_ + "'s " + strategy.title_ + ".html", "wb")

conn = httplib2.Http(".cache")

data_url = "http://eterna.cmu.edu/get_puzzles.php";

response, csv_data = conn.request(data_url, "GET")

csv_reader = csv.reader(csv_data.split("\n"))
row_num = 0


designs = []

for row in csv_reader:
    if row_num == 0 :
        tokens = row
    else :
        if(len(row) < len(tokens)):
            continue
        col_ii = 0
        design = {}
        for item in row :
            design[tokens[col_ii].lower()] = item
            col_ii += 1
        design['secstruct_elements'] = eterna_utils.get_rna_elements_from_secstruct(design['secstruct'])
        design['pairmap'] = eterna_utils.get_pairmap_from_secstruct(design['secstruct'])
			
        elements = design['secstruct_elements']
        designs.append(design)
			
    row_num += 1
    
print("Done parsing..")

for design in designs:
	args=(design,)
	ret = inv_utils.timeout(strategy.solve, args ,timeout_duration=2,default=["AAAA"])
	sequence = ret[0]
	t0=ret[1]
	if(t0 != "timeout"):
		
		natural = inv_utils.fold(sequence)[0]
		target = design['secstruct']
		if(natural == target):
			design['solve_time']=t0
		else:
			design['solve_time']="wrong solution"
	else:
		design['solve_time']=t0

## Output results
f_output.write("<strong>Title : </strong>" + strategy.title_ + "<br/>")
f_output.write("<strong>Author : </strong>" + strategy.author_ + "<br/>")

f_output.write("\n\n\n")

body_data = ""

body_data += ("<table style='border-style:none; padding-top:0px; margin:0px; width:100%'>")
body_data += ("<tr style='padding-bottom:11px'>")
body_data += ("<th>Puzzle name</th>")
body_data += ("<th>Solve Time</th>")
body_data += ("</tr>")

i_range = range(0, len(designs))
for ii in i_range: 
	body_data += ("<tr>")
	body_data += ("<td>")
	body_data += (designs[ii]['title'])
	body_data += ("</a>")
	body_data += ("</td>")
	body_data += ("<td>")
	body_data += str(designs[ii]['solve_time'])
	body_data += ("</td>")
	body_data += ("<td style='width:50px'>")
	body_data += ("</tr>")

body_data += ("</table>")

f_output.write(body_data)
f_output.close();

##don't know where to publish
#if(False and op == "PUBLISH" and strategy.publishable_):
#	headers = {'Content-type': 'application/x-www-form-urlencoded'}
#	body = {'title':strategy.title_, 'author':strategy.author_, 'url':strategy.url_, 'ktcorr':'%.3f' % final_kt_corr, 'body':body_data, 'codelen':strategy.code_length_}
#	response, data = conn.request("", "POST", headers=headers, body=urllib.urlencode(body))
#	print data

