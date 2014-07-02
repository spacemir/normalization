from nest import *
from nest import simple_plots
from matplotlib import pyplot as plt
import numpy as np
import random
import time

class HC:
	"""Hypercolumn class with only two microcolumns"""
	t_sim = 2000.0			# [ms]
	neuron_model="iaf_neuron"
	
	N_E = 30			# nr of pyramidal cells per MC
	N_I_per_mc = 5	# nr of basket cells per MC, however modelled as an independent population size N_I*n_minicolumns
	N_I = None

	defaultrate = 500.0			# if no inrates is given all MC gets this input
	default_noise_rate = 4000.0 	# distant input from 1000 neurons firing sporadically at 4.5 Hz?
	epsc_pyr_pyr = 0.0 		# [pA] amplitude exc. recurrent pyr pyr -> 0.9 mV
	epsc_pyr_bas = 16.0			# [pA] amplitude exc. from pyr to bas	-> 0.45 mV
	epsc_ext_pyr = 120.0		# [pA] amplitude exc. driving input to pyr -> 0.9 mV	
	epsc_noise_pyr = 15.0		# [pA] amplitude noise input. Gives EPSP = 0.1 mV
	epsc_noise_bas = 0.0		# [pA] amplitude noise input. Gives EPSP ?
	ipsc_bas_pyr = -340.0		# [pA] amplitude inh. postsynaptic current
	
	pyramidal_params = {"V_reset":0.0, "tau_m":20.0, "V_m":0.0, "E_L":0.0, \
				     "tau_syn":5.0, "t_ref":5.0, "C_m":1000.0, "V_th":20.0} # R_m = 20*10^6 Ohm

	basket_params = {"V_reset":0.0, "tau_m":10.0, "V_m":0.0, "E_L":0.0, \
				     "tau_syn":5.0, "t_ref":2.0, "C_m":200.0, "V_th":20.0}  # R_m = 50*10^6 Ohm
		
	d = 2.0  						# [ms] delay. same for all synapses	
		
	p_pyr_pyr = 0.40	# probability of ex.connection between pyr cells in same MC
	p_pyr_bas = 0.7			# probability of ex.connection between pyr cells and basket cells
	p_bas_pyr = 0.7			# probability of inh.connection between basket cells and pyr cells
	p_in_pyr  = 0.5			# probability of pyr cell getting exc. input
	p_in_bas = 0.3
	
	# statistics constants	
	n_cells_per_vm = 10
	n_cells_to_raster = 30
	
	debug = False

	def build(self):
		self.inputs = Create("poisson_generator",self.n_minicolumns+1)
		self.noise_input = Create("poisson_generator")						#noise input to all exc. neurons
		self.spikedetectors = Create("spike_detector", self.n_minicolumns+1)
		self.raster_spikedetector = Create("spike_detector")
		self.voltmeters = Create("voltmeter",self.n_minicolumns+1)
		
		#create populations of excitatory cells
		for n in range(self.n_minicolumns):
			mc = Create(HC.neuron_model, HC.N_E, params=HC.pyramidal_params)
			self.minicolumns.append(mc)
		
		#create population of inhibitory cells	
		self.basket_cells = Create(HC.neuron_model, HC.N_I, params=HC.basket_params)
					
		# set rate for external poisson generators
		self.set_inrates(self.inrates)
		self.built = True
		
	def set_inrates(self, inrates):	
		#set rate of external poisson input for respective minicolumn, last rate ev drive to basketcells
		for i,rate in enumerate(inrates):						#check way to avoid loop?
			SetStatus([self.inputs[i]],{"rate":rate})
		SetStatus(self.noise_input,{"rate":HC.default_noise_rate})	

	def connect(self):
		if not self.built:
			self.build()
		CopyModel("static_synapse", "pyr_excitatory", {"weight":HC.epsc_pyr_pyr, "delay":HC.d})
		CopyModel("static_synapse", "pyr_inhibitory", {"weight":HC.ipsc_bas_pyr, "delay":HC.d})
		CopyModel("static_synapse", "bas_excitatory", {"weight":HC.epsc_pyr_bas, "delay":HC.d})
		CopyModel("static_synapse", "external", {"weight":HC.epsc_ext_pyr, "delay":HC.d})
		
		#connect pyramidal cells internally and to basketcells, and basketcells to all pyramidal cells.
		#connect all ex populations to their spikedetector
		
		raster_neurons = []
		for i,mc in enumerate(self.minicolumns):
			random_1 = True
			if random_1:
				RandomConvergentConnect(mc, mc, int(HC.p_pyr_pyr*HC.N_E), model="pyr_excitatory")
				RandomConvergentConnect(mc, self.basket_cells, int(HC.p_pyr_bas*HC.N_E), model="bas_excitatory")
				RandomConvergentConnect(self.basket_cells, mc, int(HC.p_bas_pyr*HC.N_I), model="pyr_inhibitory")
				RandomDivergentConnect([self.inputs[i]], self.basket_cells, int(HC.p_in_bas*HC.N_I), model="bas_excitatory")
				RandomDivergentConnect([self.inputs[i]], mc, int(HC.p_in_pyr*HC.N_E), model="external")
			else:
				#pyr_pyr
				for pyr in mc:					
					input_neurons_pyr_pyr = random.sample(mc, int(HC.p_pyr_pyr*HC.N_E))
					ConvergentConnect(input_neurons_pyr_pyr, [pyr], model = "pyr_excitatory")
				#RandomConvergentConnect(mc, mc, int(HC.p_pyr_pyr*HC.N_E), model="pyr_excitatory")
				#input_pyr and bas
				input_neurons_pyr = random.sample(mc, int(HC.p_in_pyr*HC.N_E))
				ConvergentConnect([self.inputs[i]], input_neurons_pyr, model="external")
				input_neurons_bas = random.sample(self.basket_cells, int(HC.p_in_bas*HC.N_I))
				ConvergentConnect([self.inputs[i]], input_neurons_bas, model="bas_excitatory")					
				#pyr_bas
				for bas in self.basket_cells:
					input_neurons_pyr_bas = random.sample(mc, int(HC.p_pyr_bas*HC.N_E))
					ConvergentConnect(input_neurons_pyr_bas, [bas], model = "bas_excitatory")
				#bas_pyr
				for pyr in mc:
					input_neurons_bas_pyr = random.sample(self.basket_cells, int(HC.p_bas_pyr*HC.N_I))
					ConvergentConnect(input_neurons_bas_pyr, [pyr], model = "pyr_inhibitory")
					
			
			#Connect external noise to all neurons
			ConvergentConnect(self.noise_input, mc, HC.epsc_noise_pyr, HC.d)
			#Connect to devices
			ConvergentConnect(mc, [self.spikedetectors[i]])
			raster_neurons += mc[0:HC.n_cells_to_raster+1]
		
		raster_neurons +=self.basket_cells[0:HC.n_cells_to_raster+1]
		ConvergentConnect(raster_neurons, self.raster_spikedetector)
		
		#connect BC population to external input and spikedetector
		ConvergentConnect(self.inputs[-1:], self.basket_cells, model="bas_excitatory")
		ConvergentConnect(self.basket_cells, self.spikedetectors[-1:])
			
		#connect n cells from inh. and exc. populations to voltmeters
		for i,mc in enumerate(self.minicolumns):
			ex_cells = mc[0:HC.n_cells_per_vm]
			ConvergentConnect([self.voltmeters[i]], ex_cells)
		in_cells = self.basket_cells[0:HC.n_cells_per_vm]		
		ConvergentConnect([self.voltmeters[self.n_minicolumns]], in_cells)
			
		self.connected = True
		
	#def connect(self):
		#if not self.built:
			#self.build()
		#CopyModel("static_synapse", "pyr_excitatory", {"weight":HC.epsc_pyr_pyr, "delay":HC.d})
		#CopyModel("static_synapse", "pyr_inhibitory", {"weight":HC.ipsc_bas_pyr, "delay":HC.d})
		#CopyModel("static_synapse", "bas_excitatory", {"weight":HC.epsc_pyr_bas, "delay":HC.d})
		#CopyModel("static_synapse", "external", {"weight":HC.epsc_ext_pyr, "delay":HC.d})
		
		##connect pyramidal cells internally and to basketcells, and basketcells to all pyramidal cells.
		##connect all ex populations to their spikedetector
		#raster_neurons = []
		#print HC.N_I
		
		#for i,mc in enumerate(self.minicolumns):
			#RandomConvergentConnect(mc, mc, int(HC.p_pyr_pyr*HC.N_E), model="pyr_excitatory")
			#RandomConvergentConnect(mc, self.basket_cells, int(HC.p_pyr_bas*HC.N_E), model="bas_excitatory")
			#RandomConvergentConnect(self.basket_cells, mc, int(HC.p_bas_pyr*HC.N_I), model="pyr_inhibitory")
			
			##RandomDivergentConnect([self.inputs[i]], mc, int(HC.p_in_pyr*HC.N_E), model="excitatory")
			
			#if HC.debug:
				#print "n_inputs pyr_pyr ", int(HC.p_pyr_pyr*HC.N_E)
				#print "n_inputs_pyr_bas ", int(HC.p_pyr_bas*HC.N_E)
				#print "n_inputs bas_pyr", int(HC.p_bas_pyr*HC.N_I)
			
			#input_neurons = random.sample(mc, int(HC.p_in_pyr*HC.N_E))
			#ConvergentConnect([self.inputs[i]], input_neurons, model="external")
			##Connect noise to pyramidal cells
			#ConvergentConnect(self.noise_input, mc, HC.epsc_noise_pyr, HC.d)
			##Connect to devices
			#ConvergentConnect(mc, [self.spikedetectors[i]])
			#raster_neurons += mc[0:HC.n_cells_to_raster+1]
		
		##Connect noise to basketcells
		#ConvergentConnect(self.noise_input, self.basket_cells, HC.epsc_noise_bas, HC.d)
		
		#raster_neurons +=self.basket_cells[0:HC.n_cells_to_raster+1]
		#ConvergentConnect(raster_neurons, self.raster_spikedetector)
		
		##connect BC population to external input and spikedetector
		#ConvergentConnect(self.inputs[-1:], self.basket_cells, model="bas_excitatory")
		#ConvergentConnect(self.basket_cells, self.spikedetectors[-1:])
			
		##connect n cells from inh. and exc. populations to voltmeters
		#for i,mc in enumerate(self.minicolumns):
			#ex_cells = mc[0:HC.n_cells_per_vm]
			#ConvergentConnect([self.voltmeters[i]], ex_cells)
		#in_cells = self.basket_cells[0:HC.n_cells_per_vm]		
		#ConvergentConnect([self.voltmeters[self.n_minicolumns]], in_cells)
			
		#self.connected = True
			
	def run(self):
		if not self.connected:
			self.connect()
		Simulate(HC.t_sim)
		self.statistics()
		self.date_time = time.strftime("%d/%m/%Y") + "-" + time.strftime("%H:%M:%S")
		self.finished=True
	

	def statistics(self):
		"""Calculate exc. and inh. spiking rate for all mc and bc pop
		   and average activity for all MC:s"""
		self.spike_statistics=[]
		for i,sd in enumerate(self.spikedetectors[:-1]):
			n_spikes = GetStatus([sd], "n_events")[0]
			spike_rate = 1000.0*n_spikes/(HC.t_sim*HC.N_E)
			self.spike_statistics.append(spike_rate)
		
		self.hc_average = sum(self.spike_statistics)/float(self.n_minicolumns)
		
		#inhibitory spiking rate
		n_spikes = GetStatus([self.spikedetectors[-1]], "n_events")[0]
		spike_rate = 1000.0*n_spikes/(HC.t_sim*HC.N_I)
		self.spike_statistics.append(spike_rate)
		return self.spike_statistics

	def connection_info(self):
		text = "Incoming connections  row=from, col=to\n"
		text += "\t\tPyramidal\tBasket\tNoise\n"
		text += "Pyramidal\t%d\t\t%d\t%0.2f Hz\n" %( int(HC.N_E*HC.p_pyr_pyr), int(HC.N_E*HC.p_pyr_bas*self.n_minicolumns), HC.default_noise_rate)
		text += "Basket\t\t%d\t\t%d\t%0.2f Hz\n"	%( int(HC.N_I*HC.p_bas_pyr), 0, 0)
		return text
	
	def compact_info(self):
		info= "HYPERCOLUMN SIMULATION" + self.date_time + "\n"
		info+= "Simulation time: \%d ms\n\n" %self.t_sim
		info+= "Number of minicolumns: %d\n" % self.n_minicolumns
		info+= "N_E per minicolumn:%d\n" %HC.N_E
		info+= "N_I per minicolumn:%d\n" %HC.N_I_per_mc
		info+= "noise: " + str(self.default_noise_rate) + " ext_pyr: " + str(self.epsc_ext_pyr) + " pyr_pyr: " + \
		str(self.epsc_pyr_pyr) + " pyr_bas:" + str(self.epsc_pyr_bas) + " bas_pyr: " + str(self.ipsc_bas_pyr) + \
		"\n" + " P_in_pyr:" + str(self.p_in_pyr) + " P_pyr_pyr: " + str(self.p_pyr_pyr) + " P_pyr_bas " + \
		str(self.p_pyr_bas) + " P_bas_pyr: " + str(self.p_bas_pyr) + "p_in_bas: " + str(self.p_in_bas) +"\n"
		info+= "Average spiking rate MC:s:\n"
		for i,spikerate in enumerate(self.spike_statistics[:-1]):
			info+= "%d:  %.2f Hz (%.2f)\n" % (i, spikerate, self.inrates[i])  
		info+="Average spiking rate HC:%.2f Hz\n" % self.hc_average
		info+= "Average spiking inh. pop.%.2f Hz\n" % self.spike_statistics[-1]
		return info


	
	def simulation_info(self, all=False):	
		info = ""
		info+= "--------------------------------------\n"
		info+= "HYPERCOLUMN SIMULATION\t\t" + self.date_time + "\n"
		info+= "--------------------------------------\n"
		info+= "Simulation time: \t\t%d ms\n\n" %self.t_sim
		info+= "Number of minicolumns: \t\t%d\n" % self.n_minicolumns
		info+= "N_E per minicolumn:\t\t%d\n" %HC.N_E
		info+= "N_I per minicolumn:\t\t%d\n\n" %HC.N_I_per_mc
		info+= "noise: " + str(self.default_noise_rate) + " ext_pyr: " + str(self.epsc_ext_pyr) + " pyr_pyr: " + \
		str(self.epsc_pyr_pyr) + " pyr_bas:" + str(self.epsc_pyr_bas) +"\n" +"p_bas_pyr: " + str(self.ipsc_bas_pyr) + \
		" P_in_pyr:" + str(self.p_in_pyr) + " P_pyr_pyr: " + str(self.p_pyr_pyr) + "\n" + " P_pyr_bas " + \
		str(self.p_pyr_bas) + " P_bas_pyr: " + str(self.p_bas_pyr) + "\n" + "p_in_bas: " + str(self.p_in_bas) +"\n\n"
		info+= "\n"
		info+= "Average spiking rate MC:s:\n"
		for i,spikerate in enumerate(self.spike_statistics[:-1]):
			info+= "\t\t\t\t%d:  %.2f Hz \t\t(%.2f)\n" % (i, spikerate, self.inrates[i])  
		info+="\n"
		info+="Average spiking rate HC:\t%.2f Hz\n" % self.hc_average
		info+= "Average spiking inh. pop.\t%.2f Hz\n" % self.spike_statistics[-1]
		info+= "\n"

		if all:
			info+= "__________________________________________________\n"
			info +=self.connection_info()
			info+= "__________________________________________________\n"
		info+= "---------------------------------------\n\n"
		self.info = info
		return info	

	def voltage_traces(self, mcs=[]):
		"""Plots voltage_traces from given MC:s (or all) in each subplot, 
		   number of neurons to show for each pop. is set in HC.n_cells_to_voltmeter """
		if mcs == [] or len(mcs)>self.n_minicolumns:
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
		n_bins = int(HC.t_sim/bin_size)
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
					norm_factor = 1000.0/(bin_size*HC.N_E) 
				else: 											# ..som har annat antal neuroner
					norm_factor = 1000.0/(bin_size*HC.N_I) 
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
		

	def __init__(self, n_minicolumns = 3, inrates=None):
		self.name=self.__class__.__name__
		
		self.n_minicolumns = n_minicolumns	#number of MC:s
		if inrates == None:
			self.inrates= [HC.defaultrate]*n_minicolumns+[0.0]
		else:
			self.inrates = inrates
		
		self.minicolumns = []				#list of nodeslists for MC:s
		self.inputs = None					#list of poissongenerators for input to MC:s
		self.basket_cells = None			#basket cell population
		self.spikedetectors = None			#list of spikedetectors one for each MC and one for BCpop
		self.voltmeters = None
		self.raster_spikedetector = None
		self.date_time = None
		
		HC.N_I = HC.N_I_per_mc*self.n_minicolumns
		
		self.built = False
		self.connected = False
		self.population_names = []
		for i in range(self.n_minicolumns):
			self.population_names.append("Minicolumn " + str(i))
		self.population_names.append("Inhibitory pop")
		self.population_names = self.population_names[::-1]

		ResetKernel()	
		self.setRandomSeeds()
	
	def setRandomSeeds(self):
		n_vp = GetKernelStatus("total_num_virtual_procs")
		msd = 1004	# master seed
		msdrange1 = range(msd, msd+n_vp)	#list of seeds one for each virtual process
		pyrngs = [np.random.RandomState(s) for s in msdrange1]	#Create RNG:s
		msdrange2 = range(msd+n_vp +1, msd+1+2*n_vp)
		SetKernelStatus({'grng_seed':msd+n_vp, "rng_seeds":msdrange2}) #?
		
	def testPSPs(self):
		pass
		
