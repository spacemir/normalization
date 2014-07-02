from nest import *
from nest import simple_plots
from nest import voltage_trace
from nest import raster_plot
from matplotlib import pyplot as plt
from tests_hc32 import *
import numpy as np
import scipy.stats as stats
import random
import sys
import time

##########################################
# 3_1 + randomized neuron sizes, and synaptic weights
##########################################

class HC:
	"""Hypercolumn class with variable number of microcolumns."""
	
	#NORMTEST PARAMETERS low range
	#rates_other = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0]#, 1200.0] #input rates for other minicolumns
	#rates = np.concatenate((np.arange(0, 800, 25, "float"),np.arange(800, 4001, 250, "float")),0)
	
	# NORMTEST PARAMETERS short test low range 2
	rates_other = [200.0, 400.0, 600.0, 800.0] #, 1000.0, 1200.0] 
	rates = np.concatenate((np.arange(0, 1000, 50, "float"),np.arange(1000, 2001, 100, "float")),0)
	
	# NORMTEST MEDIUM RANGE	medium range
	#rates_other = [200.0, 1200.0, 2200.0, 3200.0] 
	#rates = np.concatenate((np.arange(0, 2000, 100, "float"),np.arange(2000, 6001, 200, "float")),0)
	
	# NORMTEST LONG RANGE	high range
	#rates_other = [200.0, 1200.0, 2200.0, 3200.0] 
	#rates = np.concatenate((np.arange(0, 2000, 100, "float"),np.arange(2000, 6001, 200, "float")),0)

	def initialize_vmem(self, nodes, c_m_mean, c_m_std):
		for gid in nodes:
			SetStatus([gid], {"C_m":np.random.normal(c_m_mean, c_m_std)})
			SetStatus([gid], {"V_m":np.random.normal(self.params["V_init_mean"], self.params["V_init_std"])})
		
	def build(self):
		self.inputs = Create("poisson_generator",self.params["n_minicolumns"]+1)
		self.pyr_noise_input = Create("poisson_generator")	#noise input to all exc neurons
		self.bas_noise_input = Create("poisson_generator")	#noise input to all inh neurons
		self.spikedetectors = Create("spike_detector", self.params["n_minicolumns"]+1)
		self.raster_spikedetector = Create("spike_detector")
		self.voltmeters = Create("voltmeter",self.params["n_minicolumns"]+1)
		
		#create populations of excitatory cells
		for n in range(self.params["n_minicolumns"]):
			mc = Create(self.params["neuron_model"], self.params["N_E"],\
				        params=self.params["pyramidal_params"])
			self.initialize_vmem(mc, self.params["C_m_mean_pyr"], self.params["C_m_std_pyr"])
			self.minicolumns.append(mc)

		#all_exc_neurons = []
		#for mc in self.minicolumns:
			#all_exc_neurons += mc
		#node_info = GetStatus(mc)
		#C_m_s = [ n["V_m"] for n in node_info]
		#plt.hist(C_m_s,20)
		#plt.show()
		
		#create population of inhibitory cells	
		self.basket_cells = Create(self.params["neuron_model"], self.N_I, \
								   params=self.params["basket_params"])
		self.initialize_vmem(mc, self.params["C_m_mean_bas"], self.params["C_m_std_bas"])
					
		# set rate for external poisson generators
		#SET IN RUN FOR NORMTEST# self.set_inrates(self.inrates) 
		self.built = True
		
	def set_inrates(self, inrates):	
		self.inrates = inrates
		#set rate of external poisson input for respective minicolumn, last rate ev drive to basketcells
		for i,rate in enumerate(inrates):						#check way to avoid loop?
			SetStatus([self.inputs[i]],{"rate":rate})
		SetStatus(self.pyr_noise_input,{"rate":self.params["pyr_noise_rate"]})
		SetStatus(self.bas_noise_input,{"rate":self.params["bas_noise_rate"]})

	def RCConnect_randomized_weights(self, source_pop, target_pop, n_sources, synapse_model, mean_weight):	
		for node in target_pop:		
			weights = np.random.normal(mean_weight, self.params["synapse_std_fraction"]*np.abs(mean_weight), n_sources)
			RandomConvergentConnect(source_pop, [node], n_sources, weight=list(weights), delay = self.params["d"], \
			model=synapse_model, options={'allow_multapses':self.params["mult"], 'allow_autapses':self.params["aut"]})
			
	def RDConnect_randomized_weights(self, source_pop, target_pop, n_targets, synapse_model, mean_weight):	
		for node in [source_pop]:		
			weights = np.random.normal(mean_weight, self.params["synapse_std_fraction"]*np.abs(mean_weight), n_targets)
			RandomDivergentConnect([node], target_pop, n_targets, weight=list(weights), delay = self.params["d"], \
			model=synapse_model, options={'allow_multapses':self.params["mult"], 'allow_autapses':self.params["aut"]})
		
	def connect(self):
		if not self.built:
			self.build()
		CopyModel("static_synapse", "pyr_external", {"weight":self.params["epsc_ext_pyr"], \
		          "delay":self.params["d"]})
		CopyModel("static_synapse", "pyr_recurrent", {"weight":self.params["epsc_pyr_pyr"], \
				  "delay":self.params["d"]}) 
		CopyModel("static_synapse", "bas_external", {"weight":self.params["epsc_ext_bas"], \
		          "delay":self.params["d"]})
		CopyModel("static_synapse", "bas_recurrent", {"weight":self.params["ipsc_bas_bas"],\
				  "delay":self.params["d"]})	
		CopyModel("static_synapse", "pyr_inhibitory", {"weight":self.params["ipsc_bas_pyr"], \
				  "delay":self.params["d"]})
		CopyModel("static_synapse", "bas_excitatory", {"weight":self.params["epsc_pyr_bas"],\
				  "delay":self.params["d"]})
		
		#connect mc:s internally and to basket pop, and basket pop to mc:s.
		#connect input to mc:s and basket pop
		raster_neurons = []
		for i,mc in enumerate(self.minicolumns):
			# driving input mc
			self.RDConnect_randomized_weights(self.inputs[i], mc, int(self.params["p_in_pyr"]*self.params["N_E"]),\
					"pyr_external", self.params["epsc_ext_pyr"])	
			# recurrent mc
			self.RCConnect_randomized_weights(mc, mc, int(self.params["p_pyr_pyr"]*self.params["N_E"]),\
					"pyr_recurrent", self.params["epsc_pyr_pyr"])
			# driving input to basket_cells as well	
			self.RDConnect_randomized_weights(self.inputs[i], self.basket_cells, int(self.params["p_in_bas"]*self.N_I),\
					"bas_external", self.params["epsc_ext_bas"])
			# mc to basket cells
			self.RCConnect_randomized_weights(mc, self.basket_cells, int(self.params["p_pyr_bas"]*self.params["N_E"]),\
					"bas_excitatory", self.params["epsc_pyr_bas"])	
			# basket cells to mc	
			self.RCConnect_randomized_weights(self.basket_cells, mc, int(self.params["p_bas_pyr"]*self.N_I),\
					"pyr_inhibitory", self.params["ipsc_bas_pyr"])
			
			# Connect external noise to mc
			ConvergentConnect(self.pyr_noise_input, mc, self.params["epsc_noise_pyr"], self.params["d"])
			# Connect to devices
			ConvergentConnect(mc, [self.spikedetectors[i]])
			raster_neurons += mc[0:self.params["n_cells_to_raster"]+1]
		
		raster_neurons +=self.basket_cells[0:self.params["n_cells_to_raster"]+1]
		ConvergentConnect(raster_neurons, self.raster_spikedetector)
		
		# bc pop to spikedetector
		ConvergentConnect(self.basket_cells, self.spikedetectors[-1:])
		# external noise to bc pop
		ConvergentConnect(self.bas_noise_input, self.basket_cells, self.params["epsc_noise_bas"], self.params["d"])
		# recurrent bc pop
		self.RCConnect_randomized_weights(self.basket_cells, self.basket_cells, int(self.params["p_bas_bas"]*self.N_I),\
					"bas_recurrent", self.params["ipsc_bas_bas"])
			
		# Connect cells from inh. and exc. populations to voltmeters
		for i,mc in enumerate(self.minicolumns):
			ex_cells = mc[0:self.params["n_cells_per_vm"]]
			ConvergentConnect([self.voltmeters[i]], ex_cells)
		in_cells = self.basket_cells[0:self.params["n_cells_per_vm"]]		
		ConvergentConnect([self.voltmeters[self.params["n_minicolumns"]]], in_cells)
			
		self.connected = True
	

		
	def run_hctest(self):
		if not self.connected:
			self.connect()
		self.set_inrates(self.inrates)
		Simulate(self.params["t_sim"])
		self.statistics()
		self.date_time = time.strftime("%d/%m/%Y") + "-" + time.strftime("%H:%M:%S")
		self.finished=True
		print "bas_bas:", int(self.params["p_bas_bas"]*self.N_I)
	
	def run_normtest(self):
		"""fI curve test. Inner for loop goes through different rate of stimulation for one mc 
		   while, input to other mc:s remain constant. Outer for loop repeats this for different
		   input rates for other mc:s"""
		   
		if not self.connected:
			self.connect()
		
		self.date_time = time.strftime("%d/%m/%Y") + "-" + time.strftime("%H:%M:%S")
		start = time.time()
		
		infos = []	
		data = {"rates_other":HC.rates_other, "rates":HC.rates,"hc_averages":[], \
				"outrates_studied_mc":[],"outrates_other_mc":[],"inh_rates":[]}
		for j, rate_other in enumerate(HC.rates_other):
			hc_averages = []
			outrates_studied_mc = []
			outrates_other_mc = []
			inh_rates = []
			for i,rate in enumerate(HC.rates):
				self.reset()
				new_rates = [rate] + (self.params["n_minicolumns"]-1)*[rate_other] + [0.0] #zero for inh pop)
				self.set_inrates(new_rates)
				#print new_rates
				Simulate(self.params["t_sim"])
				activities, hc_av = self.statistics()	
				hc_averages.append(hc_av)
				outrates_studied_mc.append(activities[0])
				outrates_other_mc.append(activities[1])
				inh_rates.append(activities[-1])
				
			data["hc_averages"].append(hc_averages)
			data["outrates_studied_mc"].append(outrates_studied_mc)
			data["outrates_other_mc"].append(outrates_other_mc)
			data["inh_rates"].append(inh_rates)
		
		end = time.time()
		print "Time for simulation: %d seconds" %(end-start)
		#print data
		return data	
			

	def reset(self):
		""" Reset devices without rebuilding network """
		#? Can this sometimes be a problem, are there other stuff I'm not resetting?
		
		for sd in self.spikedetectors:
			SetStatus([sd], {"n_events":0})
		#for mc in self.minicolumns:
		#	SetStatus(mc, {"E_m":0})
		for vm in self.voltmeters:
			pass
		self.spike_statistics = []
		self.hc_average = 0	
	
	def statistics(self):
		"""Calculates spiking rate for all mc and inh pop
		   and average activity for all MC:s"""
		self.spike_statistics=[]
		# spiking rate for mc:s
		for i,sd in enumerate(self.spikedetectors[:-1]):
			n_spikes = GetStatus([sd], "n_events")[0]
			spike_rate = 1000.0*n_spikes/(self.params["t_sim"]*self.params["N_E"])
			self.spike_statistics.append(spike_rate)

		self.hc_average = sum(self.spike_statistics)/float(self.params["n_minicolumns"])
		# inhibitory spiking rate
		n_spikes = GetStatus([self.spikedetectors[-1]], "n_events")[0]
		spike_rate = 1000.0*n_spikes/(self.params["t_sim"]*self.N_I)
		self.spike_statistics.append(spike_rate)
		
		return self.spike_statistics, self.hc_average
	
	def simulation_info(self, c_info=False):	
		info = ""
		info+= "--------------------------------------\n"
		info+= "HYPERCOLUMN SIMULATION\t\t" + self.date_time + "\n"
		info+= "--------------------------------------\n"
		info+= "Simulation time: \t\t%d ms\n\n" %self.params["t_sim"]
		info+= "Number of minicolumns: \t\t%d\n" % self.params["n_minicolumns"]
		info+= "N_E per minicolumn:\t\t%d\n" %self.params["N_E"]
		info+= "N_I per minicolumn:\t\t%d\n\n" %self.params["N_I_per_mc"]
		info+= "noise pyr, bas: " + str(self.params["pyr_noise_rate"]) + ", " + str(self.params["bas_noise_rate"]) + " ext_pyr: " + \
				str(self.params["epsc_ext_pyr"]) + " pyr_pyr: " +	str(self.params["epsc_pyr_pyr"]) + " pyr_bas:" + \
				str(self.params["epsc_pyr_bas"]) +"\n" +"p_bas_pyr: " + str(self.params["ipsc_bas_pyr"]) + " p_bas_bas: " +\
				str(self.params["p_bas_bas"]) + " p_in_pyr: " + str(self.params["p_in_pyr"]) + " p_pyr_pyr: " + str(self.params["p_pyr_pyr"]) + "\n" + " p_pyr_bas " + \
				str(self.params["p_pyr_bas"]) + " p_bas_pyr: " + str(self.params["p_bas_pyr"]) + "\n" + " p_in_bas: " +\
				str(self.params["p_in_bas"]) +"\n\n"
		info+= "\n"
		info+= "Average spiking rate MC:s:\n"
		for i,spikerate in enumerate(self.spike_statistics[:-1]):
			info+= "\t\t\t\t%d:  %.2f Hz \t\t(%.2f)\n" % (i, spikerate, self.inrates[i])  
		info+="\n"
		info+="Average spiking rate HC:\t%.2f Hz\n" % self.hc_average
		info+= "Average spiking inh. pop.\t%.2f Hz\n" % self.spike_statistics[-1]

		if c_info:
			info+= "__________________________________________________\n"
			info +=self.connection_info()
			info+= "__________________________________________________\n"
		info+= "---------------------------------------\n\n"
		self.info = info
		return info		


	def compact_simulation_info(self):
		info= "HYPERCOLUMN SIMULATION" + self.date_time + "\n"
		info+= "Simulation time: \%d ms\n\n" % self.params["t_sim"]
		info+= "Number of minicolumns: %d\n" % self.params["n_minicolumns"]
		info+= "N_E per minicolumn:\t\t%d\n" %self.params["N_E"]
		info+= "N_I per minicolumn:\t\t%d\n\n" %self.params["N_I_per_mc"]
		info+= "noise pyr, bas: " + str(self.params["pyr_noise_rate"]) + " " + str(self.params["bas_noise_rate"])  + " ext_pyr: " + str(self.params["epsc_noise_pyr"]) + " pyr_pyr: " + \
		str(self.params["epsc_pyr_pyr"]) + " pyr_bas:" + str(self.params["epsc_pyr_bas"]) + " bas_pyr: " + str(self.params["ipsc_bas_pyr"]) + \
		"\n" + " P_in_pyr:" + str(self.params["p_in_pyr"]) + " P_pyr_pyr: " + str(self.params["p_pyr_pyr"]) + " P_pyr_bas " + \
		str(self.params["p_pyr_bas"]) + " P_bas_pyr: " + str(self.params["p_bas_pyr"]) + "p_in_bas: " + str(self.params["p_in_bas"]) +"\n"
		info+= "Average spiking rate MC:s:\n"
		for i,spikerate in enumerate(self.spike_statistics[:-1]):
			info+= "%d:  %.2f Hz (%.2f)\n" % (i, spikerate, self.inrates[i])  
		info+="Average spiking rate HC:%.2f Hz\n" % self.hc_average
		info+= "Average spiking inh. pop.%.2f Hz\n" % self.spike_statistics[-1]
		return info

	def connection_info(self):
		text = "Incoming connections  row=from, col=to\n"
		text += "\t\tPyramidal\tBasket\tNoise\n"
		text += "Pyramidal\t%d\t\t%d\t%0.2f Hz\n" %( int(self.params["N_E"]*self.params["p_pyr_pyr"])), \
			int(self.params["N_E"]*self.params["p_pyr_bas"]*self.params["n_minicolumns"])
		text += "Basket\t\t%d\t\t%d\t%0.2f Hz\n"	%( int(self.N_I*self.params["p_bas_pyr"]), 0, 0)
		return text
		
	def voltage_traces(self, mcs=[]):
		"""Plots voltage_traces from given MC:s (or all) in each subplot, 
		   number of neurons to show for each pop. is set in HC.n_cells_to_voltmeter """
		if mcs == [] or len(mcs)>self.params["n_minicolumns"]:
			voltmeters = self.voltmeters
		else:
			voltmeters = []
			for i in mcs:
					voltmeters.append(self.voltmeters[i])
		subplot = len(voltmeters)*100+10
		#plt.figure()
		for i, vm in enumerate(voltmeters[::-1]):
			subplot +=1	
			ax = plt.subplot(subplot)
			plt.title(self.population_names[i])
			for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
				ax.get_xticklabels() + ax.get_yticklabels()):
				item.set_fontsize(10)
			simple_plots.voltage_trace([vm])

	def activity_histograms(self,bin_size=10):
		"""Plots spike histogram scaled to Hz for all mc:s and inh. pop"""
		#plt.figure()
		plt.title("Activity histogram")
		plt.xlabel("time [ms]")
		plt.ylabel("n_spikes")
		n_bins = int(self.params["t_sim"]/bin_size)
		n_subplots = len(self.spikedetectors)
		for i, sd in enumerate(self.spikedetectors[::-1]):
			ax = plt.subplot(n_subplots,1,i+1)
			ax.set_ylim([0, 200])
			ax.set_xlim([0, HC.t_sim])
			plt.title(self.population_names[i])
			events = GetStatus([sd],"events")[0]
			times = events["times"]
			if not len(times):
				print "0 events"
			else:
				if not i == 0: 									# Om inte inhibitory pop..
					norm_factor = 1000.0/(bin_size*self.params["N_E"]) 
				else: 											# ..som har annat antal neuroner
					norm_factor = 1000.0/(bin_size*self.N_I) 
				h_weights = [norm_factor]*len(times)
				plt.hist(times, n_bins, weights=h_weights)		# vikta staplarna for att fa i <Hz>
				plt.ylabel("[Hz]")
			self.set_font(ax)
		plt.show()

	def set_font(self,ax):
		"""Set size of font for plot with handle ax, for all text"""
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
			ax.get_xticklabels() + ax.get_yticklabels()):
			item.set_fontsize(10)


	def __init__(self, params, inrates = None):
		self.name=self.__class__.__name__
		#setKernelStatus("communicate_allgather":False)
		self.inrates = inrates
		
		self.params = params				#Load parameter dictionary from file
		
		self.minicolumns = []				#list of nodeslists. MC added here during build
		self.N_I = self.params["N_I_per_mc"]*self.params["n_minicolumns"]
		
		self.built = False
		self.connected = False
		# Save polulation names if this is needed for evaluation
		self.population_names = []
		for i in range(self.params["n_minicolumns"]):
			self.population_names.append("Minicolumn " + str(i))
		self.population_names.append("Inhibitory pop")
		self.population_names = self.population_names[::-1]
		
		ResetKernel()	
		self.setRandomSeeds()
	
	def setRandomSeeds(self):
		""" Set random seeds, supposed to work for more than one virtual process"""
		
		# If I do not explicitly stt and change master seed, this seems to be the same always
		# is that not kind of weird?
		#SetKernelStatus({"local_num_threads":4})
		n_vp = GetKernelStatus("total_num_virtual_procs")
		msd = 1004	# master seed
		msdrange1 = range(msd, msd+n_vp)	#list of seeds one for each virtual process
		pyrngs = [np.random.RandomState(s) for s in msdrange1]	#Create RNG:s
		msdrange2 = range(msd+n_vp +1, msd+1+2*n_vp)
		SetKernelStatus({'grng_seed':msd+n_vp, "rng_seeds":msdrange2}) #?
		
	def testPSPs(self):
		pass

if __name__=="__main__":
	print sys.argv[0]
	if sys.argv[1] == "n":
		normalization_test()
	elif sys.argv[1] == "h":
		hypercolumn_test()
	else:
		print "Non-existing argument"

	
	
	
		## Randomize membrane time constant for nodes
		#mu = self.params["pyramidal_params"]["tau_m"]
		#sigma = self.params["pyr_std_tau_m"]		
		#node_info = GetStatus(mc)
		#local_nodes =[(n["global_id"], n["vp"]) for n in node_info if n["local"]]
		#print local_nodes
		#for gid,vp in local_nodes: 
	 
	 
	#def run_normtest_old(self):
		#"""fI curve test. Inner for loop goes through different rate of stimulation for one mc 
		   #while, input to other mc:s remain constant. Outer for loop repeats this for different
		   #input rates for other mc:s"""
		   
		#if not self.connected:
			#self.connect()
		
		#self.date_time = time.strftime("%d/%m/%Y") + "-" + time.strftime("%H:%M:%S")
		#start = time.time()
		
		#infos = []	
		#input_output = []
		#for j, rate_other in enumerate(HC.rates_other):
			#hc_averages = []
			#inrates = []
			#outrates_studied_mc = []
			#outrates_other_mc = []
			#inh_rates = []
			#for i,rate in enumerate(HC.rates):
				#self.reset()
				#new_rates = [rate] + (self.n_minicolumns-1)*[rate_other] + [0.0] #zero for inh pop)
				#self.set_inrates(new_rates)
				##print new_rates
				#Simulate(HC.t_sim)
				#activities, hc_av = self.statistics()
				#inrates.append(rate)
				#hc_averages.append(hc_av)
				#outrates_studied_mc.append(activities[0])
				#outrates_other_mc.append(activities[1])
				#inh_rates.append(activities[-1])
			#input_output.append((rate_other, (inrates, outrates_studied_mc, \
								  #hc_averages, outrates_other_mc, inh_rates)))
			##infos.append(self.simulation_info())
		#end = time.time()
		#print "Time for simulation: %d seconds" %(end-start)
		#return input_output	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
