import sys
import svg
import re
import random

########### Plot setting #######################


BAR_WIDTH = 20
PUZZLE_MARGIN = 40
POINT_H = 10
WIDTH = 3000
HEIGHT = 1000
W_MARGIN = 100
H_MARGIN = 100
SCORE_LOWER_BOUND = 45
SCORE_STEP = 5

COLORS = {}
COLORS['Players'] = "#7777ff"
COLORS['RNAInverse'] = "#77ff77"
COLORS['NUPACK'] = "#ff7777"
COLORS['EteRNA algorithm (sparse features)'] = "#ff77ff"
COLORS['EteRNA algorithm (all features)'] = "#aa77aa"
COLORS['EteRNA algorithm (conventional)'] = "#cc77cc"
COLORS['ensemble_etc_designs'] = "#000000"

AGENTS = ['RNAInverse', 'NUPACK', 'ensemble_etc_designs', 'EteRNA algorithm (conventional)', 'EteRNA algorithm (sparse features)', 'EteRNA algorithm (all features)']
AGENTS = ['RNAInverse', 'NUPACK', 'ensemble_etc_designs', 'EteRNA algorithm (sparse features)']


BG_BAR_OPACITY = 0.6

# Set of puzzles to render
# ALL
#puzzle_nids = [17320, 24153, 26603, 321468, 382268, 405631, 437822, 444718, 452127, 465127, 479426, 497237, 502447, 507875, 525133, 535658, 633803, 669347]
# Killer plot set ("The Finger" to "The Branches")
#puzzle_nids = [17320, 24153, 26603, 321468, 382268, 405631]
# Test set (Hard Lab Run series for bots)
puzzle_nids = [479426, 497237, 502447, 507875, 525133, 669347, 741980, 809077]


########################################################






data = open("merged_out.txt",'r')
designs = []
for line in data:
	tokens = line.split("\t")
	design = {}
	design['puztitle'] = tokens[1]
	design['soltitle'] = tokens[3]
	design['score'] = float(tokens[5]) + (random.random() - 0.5)
	design['username'] = tokens[2]
	if("Round" in tokens[2]):
		design['username'] = "Player"
		round_match = re.findall("\d+",tokens[2])
		design['round'] = int(round_match[0])
	else:
		design['round'] = 1
	design['exptype'] = tokens[4]
	match = re.findall("\d+",design['puztitle'])
	design['puznid'] = int(match[0])
	designs.append(design)




last_puznid = -1
puzzles = {}
puzzle_titles = []
for nid in puzzle_nids:
	puzzle_titles.append("")


for design in designs:
	if(last_puznid != int(design['puznid'])):
		last_puznid = int(design['puznid'])
		puzzles[design['puztitle']] = {}
		puzzles[design['puztitle']]['RNAInverse'] = []
		puzzles[design['puztitle']]['NUPACK'] = []
		puzzles[design['puztitle']]['Players'] = []
		puzzles[design['puztitle']]['player_rounds'] = []
		puzzles[design['puztitle']]['EteRNA algorithm (sparse features)'] = []
		puzzles[design['puztitle']]['EteRNA algorithm (all features)'] = []
		puzzles[design['puztitle']]['EteRNA algorithm (conventional)'] = []
		puzzles[design['puztitle']]['ensemble_etc_designs'] = []
	
		for nn in range(0,len(puzzle_nids)):
			if(puzzle_nids[nn] == last_puznid):
				puzzle_titles[nn] = design['puztitle']
	
	username = (design['username'])
	soltitle = design['soltitle']
	puztitle = design['puztitle']
	
	if(float(design['score']) < 45):
		continue
	
	if(username == "InverseRNA"):
		puzzles[puztitle]['RNAInverse'].append([design['score'], design['exptype']])
	elif(username == "NUPACK"):	
		puzzles[puztitle]['NUPACK'].append([design['score'], design['exptype']])
	elif(username == "Sparse Ensemble"):	
		puzzles[puztitle]['EteRNA algorithm (sparse features)'].append([design['score'], design['exptype']])
	elif(username == "L2 ensemble"):
		puzzles[puztitle]['EteRNA algorithm (all features)'].append([design['score'], design['exptype']])
	elif(username == "Conventional ensemble"):
		puzzles[puztitle]['EteRNA algorithm (conventional)'].append([design['score'], design['exptype']])
	elif(username == "Other Ensemble"):
		puzzles[puztitle]['ensemble_etc_designs'].append([design['score'], design['exptype']])
	elif(username == "Player"):
		puzzles[puztitle]['Players'].append([design['score'], design['exptype']])
		puzzles[puztitle]['player_rounds'].append(int(design['round']))
	else:
		print "Unknown agent name %s" % username
		sys.exit(0)


for key in puzzle_titles:
	puzzle = puzzles[key]
	max_round = 0
	scores = puzzle['Players']
	rounds = puzzle['player_rounds']
	for ii in range(0,len(scores)):
		if(rounds[ii] > max_round):
			max_round = rounds[ii]
	
	puzzle['max_round'] = max_round
	
	for ii in range(1,max_round+1):
		puzzle['Players_%d' % ii] = []		

	for ii in range(0,len(scores)):
		solround = rounds[ii]
		puzzle['Players_%d' % solround].append(scores[ii])

WIDTH = W_MARGIN

for key in puzzle_titles:
	puzzle = puzzles[key]
	puzround = puzzle['max_round']
	
	WIDTH += PUZZLE_MARGIN	
	
	agents = []
	
	for ii in range(0,len(AGENTS)):
		agents.append(AGENTS[ii])
	
	
	for ii in range(1,puzround+1):
		agents.append("Players_%d" % ii)
	
	for agent in agents:
		scores = puzzle[agent]
		if len(scores) > 0:
			WIDTH += BAR_WIDTH
	WIDTH += PUZZLE_MARGIN

	
plot = svg.svg("progress_plot.svg",WIDTH + 30,HEIGHT + 30)

xaxis_y = HEIGHT - (SCORE_LOWER_BOUND) * POINT_H + H_MARGIN
xaxis_y_max = H_MARGIN - POINT_H * 2
plot.line(W_MARGIN,xaxis_y,WIDTH,xaxis_y,"#000000")
plot.line(W_MARGIN,xaxis_y_max,WIDTH,xaxis_y_max,"#000000")
plot.line(W_MARGIN,xaxis_y,W_MARGIN,xaxis_y_max,"#000000")
plot.line(WIDTH,xaxis_y,WIDTH,xaxis_y_max,"#000000")

for step in range(SCORE_LOWER_BOUND,100,SCORE_STEP):
	score = step + SCORE_STEP
	plot.line(W_MARGIN, HEIGHT - score * POINT_H + H_MARGIN, WIDTH, HEIGHT - score * POINT_H + H_MARGIN, "#949494", 1)
	plot.text(W_MARGIN - 11, HEIGHT - score * POINT_H + H_MARGIN,14,"#949494","end","%d" % score)

w_walker = W_MARGIN

agents_used = []

for key in puzzle_titles:
	
	w_walker += PUZZLE_MARGIN	
	
	puzzle = puzzles[key]
	puzround = puzzle['max_round']
	
	agents = []
	for ii in range(0,len(AGENTS)):
		agents.append(AGENTS[ii])

	for ii in range(1,puzround+1):
		agents.append("Players_%d" % ii)
	

	
	for agent in agents:
		scores = puzzle[agent]
		if("Players" in agent):
			color = COLORS['Players']
		else:
			color = COLORS[agent]
		if len(scores) > 0:
			agentname = agent
			if("Players" in agent):
				agentname = "Players"
			
			if(agentname not in agents_used):
				agents_used.append(agentname)
			
			max_score = scores[0][0]
			min_score = scores[0][0]				
			
			for score in scores:
				if(score[0] > max_score):
					max_score = score[0]
				if(score[0] < min_score):
					min_score = score[0]
			
			plot.polygon([(w_walker,HEIGHT - min_score * POINT_H + H_MARGIN), (w_walker + BAR_WIDTH, HEIGHT - min_score * POINT_H + H_MARGIN), (w_walker + BAR_WIDTH, HEIGHT - max_score * POINT_H + H_MARGIN), (w_walker, HEIGHT - max_score * POINT_H + H_MARGIN)], color, color,BG_BAR_OPACITY)
			
			radius = BAR_WIDTH/4.0
			
			for score in scores:
				if(score[1] == "NaCl"):
					plot.circle(w_walker + radius,HEIGHT - score[0] * POINT_H + H_MARGIN, radius,color,"#000000")
				else:
					plot.polygon([(w_walker + radius * 2,HEIGHT - score[0] * POINT_H + H_MARGIN - radius), (w_walker + radius * 4,HEIGHT - score[0] * POINT_H + H_MARGIN - radius) , (w_walker + radius * 4,HEIGHT - score[0] * POINT_H + H_MARGIN + radius), (w_walker + radius * 2,HEIGHT - score[0] * POINT_H + H_MARGIN + radius)], color, "#000000")
				#if(score > min_score and score < max_score):
				#	plot.line(w_walker, HEIGHT - score * POINT_H + H_MARGIN, w_walker + BAR_WIDTH, HEIGHT - score * POINT_H + H_MARGIN, "#000000", 1)
			
			
			
			w_walker += BAR_WIDTH
	w_walker += PUZZLE_MARGIN
	plot.line(w_walker, xaxis_y, w_walker, xaxis_y_max, "#949494", 1)

legends_h = (len(agents_used) + 2) * 20 + 15
		
plot.polygon([ (WIDTH - 290, xaxis_y - 20), (WIDTH - 20, xaxis_y - 20), (WIDTH - 20, xaxis_y - 20 - legends_h), (WIDTH - 290, xaxis_y - 20 - legends_h)], "#ffffff", "#000000")
for ii in range(0, len(agents_used)):
	plot.polygon([ (WIDTH - 280, xaxis_y - 10 - legends_h + 20 * ii), (WIDTH - 250, xaxis_y - 10 - legends_h + 20 * ii), (WIDTH - 250, xaxis_y -10 - legends_h + 20 * ii + 15), (WIDTH - 280, xaxis_y -10 - legends_h + 20 * ii + 15)], COLORS[agents_used[ii]], COLORS[agents_used[ii]])
	plot.text(WIDTH - 250 + 10, xaxis_y- legends_h + 2 + 20 * ii, 15,"#000000", "start", agents_used[ii])
plot.circle(WIDTH - 280 + 7.5, xaxis_y - 10 - legends_h + 20 * len(agents_used) + 7.5, 7.5,"#ffffff","#000000")
plot.text(WIDTH - 265 + 10, xaxis_y- legends_h + 2 + 20 * len(agents_used), 15,"#000000", "start", "1M NaCl")
plot.polygon([ (WIDTH - 280, xaxis_y - 10 - legends_h + 20 * (len(agents_used) + 1)), (WIDTH - 265, xaxis_y - 10 - legends_h + 20 * (len(agents_used)+1)), (WIDTH - 265, xaxis_y -10 - legends_h + 20 * (len(agents_used)+1) + 15), (WIDTH - 280, xaxis_y -10 - legends_h + 20 * (len(agents_used)+1) + 15)], "#ffffff", "#000000")
plot.text(WIDTH - 265 + 10, xaxis_y- legends_h + 2 + 20 * (len(agents_used)+1), 15,"#000000", "start", "10mM MgCl2")