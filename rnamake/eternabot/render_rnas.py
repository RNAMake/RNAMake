import sys
import svg
import re
import random
import math
import render_rna


########### Plot setting #######################


# ALL
#puzzle_nids = [17320, 24153, 26603, 321468, 382268, 405631, 437822, 444718, 452127, 465127, 479426, 497237, 502447, 507875, 525133, 535658, 633803, 669347]
# Killer plot set ("The Finger" to "The Branches")
#puzzle_nids = [17320, 24153, 26603, 321468, 382268, 405631]
# Test set (Hard Lab Run series for bots)
puzzle_nids = [479426, 497237, 502447, 507875, 525133]
#puzzle_nids = [405631]
#agents = ["Target", "Player", "InverseRNA", "NUPACK", "Conventional ensemble", "L2 ensemble", "Sparse Ensemble"] 
agents = ["Target", "Player", "InverseRNA", "NUPACK", "Sparse Ensemble"] 
#agents = ["Target", "Player", "InverseRNA", "NUPACK"]
#agents = ["Player", "InverseRNA", "NUPACK"]

NODE_R = 6
PRIMARY_SPACE = 20
PAIR_SPACE = 23

CELL_PADDING = 40
TEXT_SIZE = 40

FORCE_TARGET_SHAPE = False
USE_RIGHTWRONG_COLORS = False
RENDER_IN_LETTERS = False


# Blue - Gray - Yellow
SHAPE_COLORS = [ [0,0,255], [128, 128, 128], [255,128,0]]
# Blue - Red
RIGHTWRONG_COLORS = [ [128,128,128], [255,77,77]]

# Strong, Weak
PAIRBOND_COLOR = [[0,0,0],[255,255,255]]

##################################




list_f = open("bootstrap_data/list.list", "r")
list = list_f.read()
list_f.close()

files = list.split('\n')

puzzles = {}

for file in files:
	if len(file) > 0:
		file_f = open("bootstrap_%s" %file, 'r')
		data = file_f.read()
		file_f.close()
		
		data_lines = data.split('\n')
		puzname = data_lines[0]
		agent = data_lines[1]
		sequence = data_lines[2]
		secstruct = data_lines[3]
		target_secstruct = data_lines[4]
		shape_start = int(data_lines[5]) -1
		shape_min = float(data_lines[6])
		shape_max = float(data_lines[7])
		shape_threshold = float(data_lines[8])
		score = float(data_lines[9])
		shape = data_lines[10].split(",")
		probmat = []
		
		for ii in range(11, len(data_lines)):
			if len(data_lines[ii]) > 0:
				row = data_lines[ii].split(",")
				if len(row) != len(secstruct):
					print "wrong data row length"
					sys.exit(0)
				for jj in range(0,len(row)):
					row[jj] = float(row[jj])
				probmat.append(row)
		
		if len(probmat) != len(secstruct):
			print "wrong data matrix"
			sys.exit(0)
		
		pairmap = render_rna.get_pairmap_from_secstruct(secstruct)
		
				
		shape_obj = { "shape": shape, "shape_start":shape_start, "min": shape_min, "max": shape_max, "threshold" : shape_threshold, "score":score }
		
		
		for ii in range(0,len(shape)):
			shape[ii] = float(shape[ii])
		
		if puzname not in puzzles:
			puzzles[puzname] = {}
		
		puzzle = puzzles[puzname]
		
		if "Target" not in puzzle or FORCE_TARGET_SHAPE:
			target_shape_min = 0
			target_shape_max = 1.0
			target_shape_threshold = 0.5
			target_shape = []
			
			for ii in range(0,len(target_secstruct)):
				if(target_secstruct[ii] == "."):
					target_shape.append(1.0)
				else:
					target_shape.append(0.0)
				
				target_shape_obj = {"shape":target_shape, "shape_start" : 0, "min": target_shape_min, "max": target_shape_max, "threshold" : target_shape_threshold}
			
			target_pairs = []
			target_pairmap = render_rna.get_pairmap_from_secstruct(target_secstruct)
			for ii in range(0,len(target_pairmap)):
				if(target_pairmap[ii] > ii):
					pair_obj = {"from":ii, "to":target_pairmap[ii], "p":1.0}
					target_pairs.append(pair_obj)

			puzzle["Target"] = {"secstruct":target_secstruct, "shape":target_shape_obj, "pairs":target_pairs}
			
			if FORCE_TARGET_SHAPE:
				secstruct = target_secstruct
				pairmap = target_pairmap
					
		pairs = []
		
		for ii in range(0,len(pairmap)):
			if(pairmap[ii] > ii):
				pair_obj = {"from":ii, "to":pairmap[ii], "p":probmat[ii][pairmap[ii]]}
				pairs.append(pair_obj)


		puzzle[agent] = {"secstruct" : secstruct, "sequence" : sequence, "score":score, "shape" :shape_obj, "pairs" : pairs}
		

cell_max_w = 0
cell_max_h = 0

for puzzle_nid in puzzle_nids:
	puzname = str(puzzle_nid)
	if puzname in puzzles:
		puzzle = puzzles[puzname]
		for agent in agents:
			if agent in puzzle:
				
				renderer = render_rna.RNARenderer()
				renderer.setup_tree(puzzle[agent]["secstruct"], NODE_R,PRIMARY_SPACE, PAIR_SPACE)
				size = renderer.get_size()
				puzzle[agent]["renderer"] = renderer
				
				if(size[0] > cell_max_w):
					cell_max_w = size[0]
				if(size[0] > cell_max_h):
					cell_max_h = size[1]
				

			
cell_size = max(cell_max_w,cell_max_h)	+ CELL_PADDING * 2	

svg_w = cell_size * len(puzzle_nids)
svg_h = cell_size * len(agents)
svgobj = svg.svg("secstruct_table.svg", svg_w, svg_h)

puzzle_index = 0
for puzzle_nid in puzzle_nids:
	puzname = str(puzzle_nid)
	if puzname in puzzles:
		puzzle = puzzles[puzname]
		agent_index = 0
		for agent in agents:
			if agent in puzzle:
				colors = None
				pairs = None
				
				if 'pairs' in puzzle[agent]:
					pairs = puzzle[agent]['pairs']
					
					
					for pair in pairs:
						p = pair['p']
						color = [0,0,0]
						for ii in range(0,3):
							color[ii] = int(PAIRBOND_COLOR[0][ii] * (p + 0.2) + PAIRBOND_COLOR[1][ii] * (1-p) * 0.8)
						pair['color'] = color
					
				if 'shape' in puzzle[agent]:
					colors = []
					shape_obj  = puzzle[agent]['shape']
					shape = shape_obj['shape']
					shape_start = shape_obj['shape_start']
					shape_min = shape_obj['min']
					shape_max = shape_obj['max']
					shape_threshold = shape_obj['threshold']
					secstruct = puzzle[agent]['secstruct']
					
					for ii in range(0,len(secstruct)):
						
						color = ""
							
						if USE_RIGHTWRONG_COLORS:
							if ii <shape_start or ii >= shape_start + len(shape):
								color = RIGHTWRONG_COLORS[0]
							else:
								loop_exposed = shape[ii - shape_start] > shape_threshold * 0.25 + shape_min * 0.75
								stack_unexposed = shape[ii - shape_start] < shape_threshold
								
								if not stack_unexposed and secstruct[ii] != ".":
									color = RIGHTWRONG_COLORS[1]
								elif not loop_exposed and secstruct[ii] == ".":
									color = RIGHTWRONG_COLORS[1]
								else:
									color = RIGHTWRONG_COLORS[0]
							
						else:
							if ii < shape_start or ii >= shape_start + len(shape):
								color = SHAPE_COLORS[1]
							else:
							
								exposed = shape[ii - shape_start] > shape_threshold
							
								if exposed:
									color = SHAPE_COLORS[2]
								else:
									color = SHAPE_COLORS[0]
						
						colors.append(color)
				
				
				size = puzzle[agent]["renderer"].get_size()
				x_offset = (cell_size - size[0]) / 2
				y_offset = (cell_size - size[1]) / 2
				
				sequence = None
				if "sequence" in puzzle[agent]:
					sequence = puzzle[agent]["sequence"]
				
				puzzle[agent]["renderer"].draw(svgobj, cell_size * puzzle_index + x_offset, cell_size * agent_index + y_offset, colors, pairs, sequence, RENDER_IN_LETTERS)
				if ("score" in puzzle[agent] or agent == "Target") and TEXT_SIZE > 0:
					if (agent != "Target"):
						svgobj.text(cell_size * (puzzle_index + 1) - CELL_PADDING/2, cell_size * (agent_index) + CELL_PADDING/2 + TEXT_SIZE/2, TEXT_SIZE, "#000000", "end", str(puzzle[agent]['score']))
					if puzzle_index == 0 :
						svgobj.text(cell_size * puzzle_index + CELL_PADDING/2, cell_size * agent_index + CELL_PADDING/2 + TEXT_SIZE/2, TEXT_SIZE, "#000000", "left", agent)
			agent_index += 1		
	puzzle_index += 1				
	

for ii in range(0,len(agents) + 1):
	svgobj.line(0,ii*cell_size,svg_w,ii*cell_size,"#555555")

for ii in range(0,len(puzzle_nids)+1):
	svgobj.line(ii*cell_size,0,ii*cell_size,svg_h,"#555555")
	
print "SVG size (%d, %d)" % (cell_size * len(puzzle_nids), cell_size * len(agents))	