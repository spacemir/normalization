from nest import *
import time
from collections import Counter
from HC_ModelA import HC
from HC_STD import HC_STD
import numpy as np


def hc_test(hypercolumn, inrates, bas_inrate = None):
	""" Simulates the network for t_sim ms with the specified inrates"""
	if not hypercolumn.connected:
		hypercolumn.connect()
	hypercolumn.set_inrates(inrates, bas_inrate)
	hypercolumn.simulate(hypercolumn.params["t_sim"])
	statistics(hypercolumn)
	hypercolumn.date_time = time.strftime("%d/%m/%Y") + "-" + time.strftime("%H:%M:%S")
	hypercolumn.finished = True
	
def io_test(hypercolumn):
	"""IO curve test. Inner for loop goes through different rate of stimulation for one mc 
	   while, input to other mc:s remain constant. Outer for loop repeats this for different
	   input rates for other mc:s. All simulation data is saved in the dictionary "data"
	   and returned to the caller."""   
	if not hypercolumn.connected:
		hypercolumn.connect()	
	
	hypercolumn.date_time = time.strftime("%d/%m/%Y") + "-" + time.strftime("%H:%M:%S")
	start = time.time()
	
	# dictionary for saving simulation data
	data = {"rates_other":hypercolumn.rates_other, "rates":hypercolumn.rates,"hc_averages":[], \
			"outrates_studied_mc":[],"outrates_other_mc":[],"inh_rates":[]}
	
	# perform simulations
	for j, rate_other in enumerate(hypercolumn.rates_other):
		hc_averages = []
		outrates_studied_mc, outrates_other_mc, inh_rates = [], [], []
		# rates for other minicolumns, c_i = percentage of their average input for each mc
		other_total = rate_other * hypercolumn.params["n_minicolumns"]
		c_is = np.linspace(0.2, 0.8, hypercolumn.params["n_minicolumns"]-1)\
				/sum(np.linspace(0.2,0.8,hypercolumn.params["n_minicolumns"]))
		rates_other = c_is*other_total
		print rates_other
		
		# run simulation for different inrates to studied mc
		for i,rate in enumerate(hypercolumn.rates):
			hypercolumn.reset()
			new_rates = np.concatenate((np.array([rate]),rates_other))
			hypercolumn.set_inrates(new_rates)
			hypercolumn.simulate(hypercolumn.params["t_sim"])
			activities, hc_av = statistics(hypercolumn)	
			hc_averages.append(hc_av)
			outrates_studied_mc.append(activities[0])
			outrates_other_mc.append(activities[2])
			inh_rates.append(activities[-1])
			
		data["hc_averages"].append(hc_averages)
		data["outrates_studied_mc"].append(outrates_studied_mc)
		data["outrates_other_mc"].append(outrates_other_mc)
		data["inh_rates"].append(inh_rates)
	
	end = time.time()
	print "Time for IO test simulation: %d seconds" %(end-start)
	hypercolumn.finished = True
	return data	

def statistics(hypercolumn):
	"""Calculates average spike rate for all mc:s and inh pop
	   and average activity in hypercolumn"""
	hypercolumn.spike_statistics=[]
	t_d = 50.0
	# spiking rate for mc:s
	for i,sd in enumerate(hypercolumn.spikedetectors[:-1]):
		time = GetKernelStatus("time")		
		n_spikes = getEvents([sd], time-(hypercolumn.params["t_sim"]-t_d), time)
		spike_rate = 1000.0*n_spikes/((hypercolumn.params["t_sim"]-t_d)*hypercolumn.params["N_E"])
		hypercolumn.spike_statistics.append(spike_rate)

	hypercolumn.hc_average = sum(hypercolumn.spike_statistics)/float(hypercolumn.params["n_minicolumns"])
	# inhibitory spiking rate
	n_spikes = GetStatus([hypercolumn.spikedetectors[-1]], "n_events")[0]
	spike_rate = 1000.0*n_spikes/((hypercolumn.params["t_sim"]-t_d)*hypercolumn.params["N_I"])
	hypercolumn.spike_statistics.append(spike_rate)
	
	return hypercolumn.spike_statistics, hypercolumn.hc_average

def getEvents(spikedetector, t_start, t_stop):
	""" Extract events within timeperiod ]t_start,t_stop]"""
	events = GetStatus(spikedetector, "events")[0]["times"]
	events_stripped = [x for x in events if x>t_start and x <= t_stop]
	return len(events_stripped)

def spike_frequency_analysis(hypercolumn, neuron_pop_nr = 0):
	"""Returns spike frequencies and neuron id:s for single neurons connected to spikedetector"""	
	k = neuron_pop_nr
	spikedetector = hypercolumn.spikedetectors[k]
	if k == hypercolumn.params["n_minicolumns"]:
		neuron_pop = hypercolumn.basket_cells
	else:
		neuron_pop = hypercolumn.minicolumns[k]
		
	senders = GetStatus([spikedetector], "events")[0]["senders"]
	counter = Counter(senders)
	keys = counter.keys()
	spikefrequencies = [0]*len(neuron_pop)

	for i,n in enumerate(neuron_pop):
		spikefrequencies[i] = counter[n]*1000.0/hypercolumn.params["t_sim"]
	
	return spikefrequencies, neuron_pop		

def neuron_param_histogram(prop, neuron_pop):
	"""plots a histogram of the property prop"""
	node_info = GetStatus(neuron_pop)
	data = [ n["V_m"] for n in node_info]
	plt.hist(data,20)
	plt.show()

def simulation_info(hypercolumn):	
	info = ""
	info+= "--------------------------------------\n"
	info+= "HYPERCOLUMN SIMULATION\t\t" + hypercolumn.date_time + "\n"
	info+= "--------------------------------------\n"
	info+= "Simulation time: \t\t%d ms\n\n" %hypercolumn.params["t_sim"]
	info+= "Number of active minicolumns: \t\t%d\n" % hypercolumn.params["n_minicolumns"]
	info+= "N_E per minicolumn:\t\t%d\n" %hypercolumn.params["N_E"]
	info+= "N_I total:\t\t%d\n\n" %hypercolumn.params["N_I"]
	info+= "noise pyr, bas: " + str(hypercolumn.params["pyr_noise_rate"]) + ", " + str(hypercolumn.params["bas_noise_rate"]) + " ext_pyr: " + \
			str(hypercolumn.params["epsp_ext_pyr"]) + " pyr_pyr: " +	str(hypercolumn.params["epsp_pyr_pyr"]) + " pyr_bas:" + \
			str(hypercolumn.params["epsp_pyr_bas"]) +"\n" +"p_bas_pyr: " + str(hypercolumn.params["ipsp_bas_pyr"]) + " p_bas_bas: " +\
			str(hypercolumn.params["p_bas_bas"]) + " p_in_pyr: " + str(hypercolumn.params["p_in_pyr"]) + " p_pyr_pyr: " + str(hypercolumn.params["p_pyr_pyr"]) + "\n" + " p_pyr_bas " + \
			str(hypercolumn.params["p_pyr_bas"]) + " p_bas_pyr: " + str(hypercolumn.params["p_bas_pyr"]) + "\n" + " p_in_bas: " +\
			str(hypercolumn.params["p_in_bas"]) +"\n\n"
	info+= "\n"
	if isinstance(hypercolumn,HC_STD):
		inrate_scaling = hypercolumn.params["N_E_4"]/2.0
	else:
		inrate_scaling = 1.0
	info+= "Average spiking rate MC:s:\n"
	for i,spikerate in enumerate(hypercolumn.spike_statistics[:-1]):
		info+= "\t\t\t\t%d:  %.2f Hz \t\t(%.2f)\n" % (i, spikerate, hypercolumn.inrates[i]*inrate_scaling)  
	info+="\n"
	info+="Average spiking rate HC:\t%.2f Hz\n" % hypercolumn.hc_average
	info+= "Average spiking inh. pop.\t%.2f Hz\n" % hypercolumn.spike_statistics[-1]
	hypercolumn.info = info
	return info		

