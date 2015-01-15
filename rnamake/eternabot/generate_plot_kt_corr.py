import csv
import httplib2
import sys
import re
import urllib
import random
import eterna_utils
import test_utils
import os
from scipy.optimize import fmin
import numpy
import scipy.linalg




def generate_plot(designs, strategies, ytics, ylabels, name):

	plotout = open("strategy_plot.dat", "w+")	
	strategy_walker = 0
	
	xtics = []		
	
	for strategy in strategies:
		
		xtics.append(strategy.title_)
		
		for ii in range(0, len(designs)):
			design = designs[ii]
			kt_score = 0
			kt_found = False
			for jj in range(1,len(ytics)):
				if(ytics[jj] > ii):
					 kt_score = strategy.kt_local_scores[jj-1]
					 kt_found = True
					 break
			if(kt_found == False):
				kt_score = strategy.kt_local_scores[len(ytics)-1]
			plotout.write("%d %d %f\n" % (strategy_walker, ii, kt_score))
		
		plotout.write("\n");	
	
		for ii in range(0, len(designs)):
			design = designs[ii]
			kt_score = 0
			kt_found = False
			for jj in range(1,len(ytics)):
				if(ytics[jj] > ii):
					 kt_score = strategy.kt_local_scores[jj-1]
					 kt_found = True
					 break
			if(kt_found == False):
				kt_score = strategy.kt_local_scores[len(ytics)-1]
			plotout.write("%d %d %f\n" % (strategy_walker+1, ii, kt_score))
			
		plotout.write("\n");		
		strategy_walker+=1		
	plotout.close()

	dis_name = "%s_plot_discrete.dem" % name
	con_name = "%s_plot_continuous.dem" % name

	plotdem_dis = open(dis_name,"w+")
	plotdem_con = open(con_name, "w+")
	
	plotdem_dis.write("set terminal png size 1600, 1600\n")
	plotdem_dis.write("set pm3d map\n")
	plotdem_con.write("set terminal png size 1600, 1600\n")
	plotdem_con.write("set pm3d map\n")
	
	if(len(xtics) >0):
		plotdem_dis.write("set format x \"\"\n")
		plotdem_dis.write("set xtics 1\n")
		plotdem_con.write("set format x \"\"\n")
		plotdem_con.write("set xtics 1\n")
		for ii in range(0,len(xtics)):
			plotdem_dis.write("set label \"%s\" font 'arial, 9' rotate by -45 at %f,-2\n" % (xtics[ii],(ii+0.5)))
			plotdem_con.write("set label \"%s\" font 'arial, 9' rotate by -45 at %f,-2\n" % (xtics[ii],(ii+0.5)))
	if(len(ytics) == 0):
		plotdem_dis.write("unset ytics\n")
		plotdem_con.write("unset ytics\n")
	else:
		plotdem_dis.write("set ytics font 'arial,9'\n")
		plotdem_dis.write("set ytics (")
		plotdem_con.write("set ytics font 'arial,9'\n")
		plotdem_con.write("set ytics (")
		for ii in range(0,len(ytics)):
			if ii>0:
				plotdem_dis.write(", ")
				plotdem_con.write(", ")
			plotdem_dis.write("\"\" %d" % ytics[ii])
			plotdem_con.write("\"\" %d" % ytics[ii])
		plotdem_dis.write(")\n")
		plotdem_con.write(")\n")	
			
		for ii in range(0,len(ytics)):
			ytic = 0
			if(ii < len(ytics)-1):
				ytic = (ytics[ii] + ytics[ii+1]) /2.0
			else:
				ytic = (ytics[ii] + len(designs)) /2.0
			plotdem_dis.write("set label \"%s\" right font 'arial,9' at -0.5,%f\n" % (ylabels[ii],ytic))
			plotdem_con.write("set label \"%s\" right font 'arial,9' at -0.5,%f\n" % (ylabels[ii],ytic))
			
	plotdem_dis.write("set grid lt 1 lc rgb '#000000'\n")
	plotdem_dis.write("set grid xtics\n")
	
	plotdem_con.write("set grid lt 1 lc rgb '#000000'\n")
	plotdem_con.write("set grid xtics\n")

	plotdem_dis.write('set palette model RGB defined (0 "black", 1 "blue", 2 "green", 3 "orange", 4 "red")\n')
	plotdem_con.write('set palette model RGB defined (0 "black", 1 "blue", 2 "green", 3 "orange", 4 "red")\n')
	
	plotdem_dis.write("splot \"strategy_plot.dat\" notitle\n")
	plotdem_con.write("splot \"strategy_plot.dat\" notitle\n")
	
	plotdem_dis.close()
	plotdem_con.close()
	
	os.system("C:/gnuplot440/binary/gnuplot.exe %s > %s_out_discrete.png" % (dis_name, name))
	os.system("C:/gnuplot440/binary/gnuplot.exe %s > %s_out_continuous.png" % (con_name, name))




print "QUERYING server"
designs = eterna_utils.get_synthesized_designs_from_eterna_server()
print "DONE parsing"

designs = sorted(designs, key=lambda design: (design['puznid'], design['score']), reverse=False)


path="strategies/"
dirList=os.listdir(path)
strategies = []


ytics = []
ylabels = []
last_puznid = 0
for ii in range(0, len(designs)):
	design = designs[ii]
	if(last_puznid != design['puznid']):
		last_puznid = design['puznid']
		ytic = ii
		#if(ytic <0):
			#ytic = 0
		ytics.append(ytic)
		ylabels.append(design['puztitle'])


for fname in dirList:
	if(re.search('\.py$',fname) and re.search('ensemble',fname) == None and re.search('nupack',fname) == None):
		strategy = eterna_utils.load_strategy_from_file('strategies\\' + fname)
		strategy.kt_score = test_utils.calculate_kt_corr(strategy.default_params_,designs,strategy)
		strategy.kt_local_scores = []
		for ii in range(0, len(ytics)):
			if ii < len(ytics) -1:
				local_designs = designs[ytics[ii]:ytics[ii+1]]
			else:
				local_designs = designs[ytics[ii]:len(designs)]

			strategy.kt_local_scores.append(test_utils.calculate_kt_corr(strategy.default_params_,local_designs,strategy))
		strategies.append(strategy)


strategies = sorted(strategies, key=lambda strategy: strategy.kt_score, reverse=True)
generate_plot(designs,strategies,ytics,ylabels,"kt_puzzlesort")		

designs = sorted(designs, key=lambda design: (design['puznid'], design['score']), reverse=False)

for ii in range(0,len(ytics)):
	strategies = sorted(strategies, key=lambda strategy: strategy.kt_local_scores[ii], reverse=False)
	title = ylabels[ii]
	title = title.replace(' ','_')
	title = title.replace(':','')
	print "Graphing %s" % title
	generate_plot(designs,strategies,ytics,ylabels,"kt_puzzlesort_%s" % title)
