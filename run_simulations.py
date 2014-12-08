from nest import *
from matplotlib import pyplot as plt
from scipy.optimize import leastsq
import numpy as np
import random
import sys
import time
import hc_visualization as vis
import utility_functions as uf
import data_analysis as da
import mutual_entropy as me
import hc_utils as hcu


#########################################################################
# Runs the test specifiec by user as command line argument
# h = hypercolumn_test (single input vector)
# io = IO test
# fir = FIRtest
# 2mc = test studying output from a network with only two activ minicolumns
#	    as a function of the inputs to those minicolumns.
#########################################################################


#from HC_ModelA import HC							# Hypercolumn to test
from HC_STD import HC_STD as HC

#from parameters_A_standard import params			# Choose parameter set, loaded as parameter dictionary
#from parameters_A_tuned_variability import params
from parameters_STD_tuned import params

Hyperkolumn = HC	

fir_test_params = {"1234":[0.1, 0.2, 0.3, 0.4], "1200":[1.0,2.0,0,0], "1111":[1.0,1.0,1.0,1.0], "average":[0.1, 0.2, 0.3, 0.4]}

def hypercolumn_test():
	
	inputs = [200.0, 400.0, 600.0, 800.0]

	hypercolumn = Hyperkolumn(params)
	hypercolumn.build(), hypercolumn.connect()
	hcu.hc_test(hypercolumn, inputs)
	view = vis.Visualization(hypercolumn)
	print hcu.simulation_info(hypercolumn)

	# String represents what plots to see h = histograms, v = voltages traces, 
	# s = spike frequency diagram, r = raster plot
	view.visualize("hsrv")

def IO_test(file_only = False, datainpic = False):
	io, fit = True, False
	file_name = "io_test"

	if not file_only:
		hypercolumn = HC(params)
		hypercolumn.build(), hypercolumn.connect()
		data = hcu.io_test(hypercolumn)
		date = hypercolumn.date_time
		data["date"] = date
		data["sim_info"] = hcu.simulation_info(hypercolumn)
		uf.save_obj(data, file_name)
	else:
		data = uf.load_obj(file_name)

	# Plots IO curves 
	fig, ax = plt.subplots(1)
	rates_in = data["rates"]
	colors = ['b','r','g','m','c','b','r','g','m','c','b','r','g']
	for i in range(len(data["rates_other"])):
		plt.plot(rates_in, data["outrates_studied_mc"][i], color=colors[i])			# output aktuell mc
		plt.plot(rates_in, data["hc_averages"][i], "--", color=colors[i])			# medel hc
		plt.plot(rates_in, data["outrates_other_mc"][i],":", color=colors[i])	    # output annan mc
	plt.ylabel("Output (Hz)")
	plt.xlabel("Input studied mc (Hz)")
	plt.show()

def FIR_test(file_only = False, histogram = True):
	n_mc, n_trials = 4,1	# n_mc = number of minicolumn, n_trials = number of simulations to average output over
	file_name = "fir_test"
	test= "1234"
	c = np.array(fir_test_params[test])
	c = c/float(sum(c))

	# create set of input vectors
	low, high, inc = 50, 10000, 1.3
	x = c*low*n_mc
	all_inputs = [x]
	while sum(x)/len(x)<high:
		x=x*inc
		all_inputs.append(x)
		
	if not file_only:
		hypercolumn = Hyperkolumn(params)
		hypercolumn.build(), hypercolumn.connect()
		results = np.zeros((len(all_inputs),len(all_inputs[0])+1))

		for k in range(n_trials):
			statistics = np.zeros( (len(all_inputs),len(all_inputs[0])+1) )

			for i,inputs in enumerate(all_inputs):
				hcu.hc_test(hypercolumn, inputs)
				hcu.statistics(hypercolumn)
				statistics[i] = hypercolumn.spike_statistics
				hypercolumn.reset()

			results = results + np.array(statistics)
			print results

		results = results/n_trials
		data = (all_inputs, results)
		print "Data", data
		uf.save_obj(data, file_name)

	else:
		data = uf.load_obj(file_name)
		print "Data:", data
		all_inputs, results = data[0], data[1]

	for i in range(len(all_inputs)):
		print results[i], np.array(all_inputs[i])

	# Best sequence that pass IO-test
	best_seq, start, stop, best_q = uf.FIR_score_and_graph(np.array(all_inputs),results, \
	                                                 2, test, 250, histogram)

	best_seq, start, stop, best_q = uf.FIR_score_and_graph(np.array(all_inputs),results, \
	                                                 2, "average", 250, False)
	print "Best average", best_qa, best_seqa
	return best_q, best_seq, best_qa

def hypercolumn_2mc_test(from_file = True):
	min_x, max_x, n_x = (000, 5000, 10)
	file_name = "2mc_test"
	if not from_file:
		X1 = np.linspace(min_x,max_x,n_x)
		X2 = np.linspace(min_x,max_x,n_x)
		hypercolumn = Hyperkolumn(params)
		hypercolumn.build(), hypercolumn.connect()

		Y1 = np.zeros((len(X1),len(X2)))
		Y2 = np.zeros((len(X1),len(X2)))
		for row,x1 in enumerate(X1):
			for col,x2 in enumerate(X2):

				hcu.hc_test(hypercolumn,[x1,x2,0,0])
				Y1[row][col] = hypercolumn.spike_statistics[0]
				Y2[row][col] = hypercolumn.spike_statistics[1]
				print "I = (", x1,",", x2, ")"
				print "Spikestatistics", hypercolumn.spike_statistics
				hypercolumn.reset()

		results = (X1, X2, Y1, Y2)
		uf.save_obj(results, file_name)
		print "X1:", X1
		print "X2:", X2
	else:
		results = uf.load_obj(file_name)

	illustration_2D(results)

def illustration_2D(results):
	X1,X2,Y1,Y2 = results		
#	uf.surface_3d(X1,X2,Y1, 200)
#	uf.surface_3d(X1, X2, Y2, 200)
	uf.surface_3d(X1, X2, (Y1+Y2)/4,100)
	R =  uf.divide_without_infinity(Y1,Y2)
#	uf.surface_3d(X1, X2, R, 10)
	uf.contour_2d(X1, X2, R)


if __name__=="__main__":
	np.set_printoptions(precision=2), np.set_printoptions(suppress=True)
	from_file = False
	if len(sys.argv) > 2 and sys.argv[2] == "file":
			from_file = True
	if sys.argv[1] == "io":
		IO_test(from_file)
	elif sys.argv[1] == "h":
		hypercolumn_test()
	elif sys.argv[1] == "fir":
		FIR_test(from_file)
	elif sys.argv[1] == "2mc":
		hypercolumn_2mc_test(from_file)
	else:
		raise RuntimeError("Non-existing argument, use h/io/fir/2mc")

