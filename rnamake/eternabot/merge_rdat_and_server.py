import csv
import httplib2
import sys
import random
import math
import vienna_parameters
import os
import re
import imp
import eterna_utils


def find_from_rdat(nid, rdats):
	results = None
	for rdat in rdats:
		if rdat['ID'] == nid:
			if results == None:
				results = []
			results.append(rdat)

	return results

def find_from_server(nid, designs):
	for design in designs:
		if design['solnid'] == nid:
			return design
	
	return None



dirList=os.listdir("rdats")

rdat_designs = []

for fname in dirList:
	
	if re.search("rdat", fname) == None:
		continue
	
	filename = "rdats/%s" % fname
	f = open(filename,'r')
	
	annot_index = 0
	
	for line in f:
		if re.search("ANNOTATION_DATA", line) == None:
			continue
		
		annot_index += 1
		
		design = {"from" : filename, "annot_index" : annot_index }
		
		match = re.findall(r'([^\s:]+):([^\s:]+):([^\s:]+)',line)
		if match != None:
			for tuple in match:
				design[tuple[1]] = tuple[2]
				
		
		match = re.findall(r'([^\s:]+):([^\s:]+):([^\s:]+):([^\s:]+)',line)	
		if match != None:
			for tuple in match:
				design[tuple[2]] = tuple[3]

		failure = False

		if 'design_name' not in design:
			print "design_name not in. Unzi!"
			sys.exit(0)
					
		if 'ID' not in design:
			print "ID not in. Unzi!"
			sys.exit(0)
		else:
			design['ID'] = int(design['ID'])
			
		if 'subround' not in design:
			print "subround not in. Unzi!"
			sys.exit(0)
		
		if 'target' not in design:
			design['target'] = "Unknown puzzle"
		
		if 'EteRNA_score' not in design:
			print "EteRNA_score not in. Unzi!"
			sys.exit(0)
		else:
			design['EteRNA_score'] = float(design['EteRNA_score'])
				
		if 'MgCl2' not in design and 'NaCl' not in design and "FMN" not in design:
			print "chemical not in. Unzi!" 
			sys.exit(0)
		
		if 'min_SHAPE' not in design or 'max_SHAPE' not in design or 'threshold_SHAPE' not in design:
			print "min, max or threshold not in. Unzi!"
			failure = True
		else:
			design['min_SHAPE'] = float(design['min_SHAPE'])
			design['max_SHAPE'] = float(design['max_SHAPE'])
			design['threshold_SHAPE'] = float(design['threshold_SHAPE'])
		
		if "FMN" not in design and not failure:
			rdat_designs.append(design)

	f.close()	
	
print "Total : %d" % (len(rdat_designs))

temp_designs = eterna_utils.get_synthesized_designs_from_eterna_server(False,"http://eterna.cmu.edu/eterna_get_synthesized.php?all=1",True)
designs = []
for design in temp_designs:
	if int(design['score']) > 20:
		designs.append(design)


print "Total (server) %d : " % len(designs)

print "Missing in RDAT"
for design in designs:
	nid = int(design['solnid'])
	design['solnid'] = nid
	score = float(design['score'])
	design['score'] = score
	
	if find_from_rdat(nid, rdat_designs) == None and score > 20:
		print "Missing %d %s" % (design['solnid'], design['puztitle'])
			
print "\n\n"
print "Missing in server"
for design in rdat_designs:
	if find_from_server(design['ID'], designs) == None:
		design['ignore'] = True
		print "Missing %d %s" % (design['ID'], design['target'])

print "\n\n"		
print "More than 2 scores"
for design in designs:
	results =  find_from_rdat(design['solnid'], rdat_designs)
	if results != None and len(results) > 0:
		
		froms = []

		for result in results:			
			if result['from'] not in froms:
				froms.append(result['from'])
		
		from_to_count = {}
		
		for fr in froms:
			num_f_count = 0
			for result in results:
				if result['from'] == fr:
					num_f_count += 1
			from_to_count[fr] = num_f_count
			
		
		score_diff_smallest = -1
		score_diff_signed = -1
		score_to_match  = design['score']
		closest_from = None				
		for result in results:	
			diff = abs(result['EteRNA_score'] - score_to_match)
			if diff < score_diff_smallest or score_diff_smallest < 0:
				score_diff_smallest = diff
				score_diff_signed = result['EteRNA_score'] - score_to_match
				closest_from = result['from']
		
		if score_diff_smallest > 5:
			print "C8 Unzi %f %f %s %s" % (score_to_match, score_diff_signed, design['solnid'], design['puztitle'])

		if len(results) %2 != 0:
			print "SS8 Unzi %d %s %s" % (len(results), design['solnid'], design['puztitle'])
		
		if len(results) > 3 and len(froms) < 2:
			print "Cyang Unzi!"
			
		best_from_count = 0
		for result in results:
			if result['from'] == closest_from:
				best_from_count += 1
		
		if best_from_count != 2:
			print "Anwa unzi! %d" % best_from_count
		
		for result in results:
			if result['from'] != closest_from:
				result['ignore'] = True
			


# Final Set
final_rdats = []
for rdat in rdat_designs:
	if 'ignore' in rdat:
		continue
	final_rdats.append(rdat)

# Must have N (from server) * 2 -4 (missing 2 designs) - 2 (2 designs with only 1 data points) data points 	
print "Final total %d" % len(final_rdats)

final_f = open("merged_out.txt","w+")
final_print = 0
design_index = 1
for design in designs:
	uid = int(design['uid'])
	name = ""
	if uid == 24553:
		name = "NUPACK"
	elif uid == 24195:
		name = "InverseRNA"
	elif uid == 26533:
		soltitle = design['soltitle'].lower()
		if re.search("sparse", soltitle) != None:
			name = "Sparse Ensemble"
		elif re.search("l2", soltitle) != None:
			name = "L2 ensemble"
		elif re.search("conventional", soltitle) != None:
			name = "Conventional ensemble"
		else:
			name = "Other Ensemble"
	else:
		name = "Round %s" % design['round'] 	
	rdats = find_from_rdat(design['solnid'], final_rdats)
	if rdats != None:
		for rdat in rdats:
			chem = ""
			if "NaCl" in rdat:
				chem = "NaCl"
			elif "MgCl2" in rdat:
				chem = "MgCl2"
			else:
				print "Fuck Unzi!"
				sys.exit(0)
			final_f.write("%d\t%s\t%s\t%s\t%s\t%f\t%s\t%d\t%f\t%f\t%f\n" % (design_index,design['puznid'], name, design['soltitle'], chem, rdat['EteRNA_score'], rdat['from'], rdat['annot_index'], rdat['min_SHAPE'], rdat['max_SHAPE'], rdat['threshold_SHAPE']))		
		final_print += len(rdats)
		design_index += 1
final_f.close()

puzzle_bests = {}
for design in designs:
	uid = int(design['uid'])
	name = ""
	if uid == 24553:
		name = "NUPACK"
	elif uid == 24195:
		name = "InverseRNA"
	elif uid == 26533:
		soltitle = design['soltitle'].lower()
		if re.search("sparse", soltitle) != None:
			name = "Sparse Ensemble"
		elif re.search("l2", soltitle) != None:
			name = "L2 ensemble"
		elif re.search("conventional", soltitle) != None:
			name = "Conventional ensemble"
		else:
			name = "Other Ensemble"
	else:
		name = "Player"
		
	rdats = find_from_rdat(design['solnid'], final_rdats)
	if rdats != None:
		for rdat in rdats:
			rdat_string = "%d\t%s\t%s\t%s\t%s\t%f\t%s\t%d\t%f\t%f\t%f\n" % (design_index,design['puznid'], name, design['soltitle'], chem, rdat['EteRNA_score'], rdat['from'], rdat['annot_index'], rdat['min_SHAPE'], rdat['max_SHAPE'], rdat['threshold_SHAPE'])
			if design['puznid'] not in  puzzle_bests:
				puzzle_bests[design['puznid']] = {}
						
			puzzle = puzzle_bests[design['puznid']]

			if name in puzzle:
				skip = design['puznid'] == "405631" and design['soltitle'] == "Cyborg 1"
				
				if puzzle[name]['score'] < rdat['EteRNA_score'] and not skip:
					print "%f %f" % (puzzle[name]['score'], rdat['EteRNA_score'])
					puzzle[name]['score'] = rdat['EteRNA_score']
					puzzle[name]['rdat_string'] = rdat_string
			else:
				puzzle[name] = {}
				puzzle[name]['score'] = rdat['EteRNA_score']
				puzzle[name]['rdat_string'] = rdat_string
				
			
			
final_mf = open("merged_out_matlab.txt", "w+")
final_mf.write("index\tpuznid\tuser\tsoltitle\tchemical\tscore\tfrom\tfromindex\tSHAPEmin\tSHAPEmax\tSHAPEthreshold\n")
for puznid in puzzle_bests:
	puzzle = puzzle_bests[puznid]
	for name in puzzle:
		agent = puzzle[name]
		final_mf.write(agent['rdat_string'])
final_mf.close()


print "Final lines %d" % final_print