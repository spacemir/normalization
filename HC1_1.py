from nest import *
from nest import simple_plots
from matplotlib import pyplot as plt
import time

class HC:
	"""Class for collecting nodes and data about hypercolumn"""
	t_sim = 500.0			# [ms]
	neuron_model="iaf_neuron"
	
	N_E = 30			# nr of pyramidal cells per MC
	N_I = 4 			# nr of basket cells per MC, however modelled as an independent population size N_I*n_minicolumns

	defaultrate = 5000.0
	epsc_recurrent = 50.0 		# [pA] amplitude exc. postsynaptic current 
	ipsc_recurrent = -150.0		# [pA] amplitude inh. postsynaptic current
	
	basket_params = {"V_reset":0.0, "tau_m":5.0, "V_m":0.0, "E_L":0.0, \
				     "tau_syn":2.0, "t_ref":2.0, "C_m":200.0, "V_th":20.0}  # R_m = 40*10^6 Ohm
				    
	pyramidal_params = {"V_reset":0.0, "tau_m":20.0, "V_m":0.0, "E_L":0.0, \
				     "tau_syn":2.0, "t_ref":5.0, "C_m":1000.0, "V_th":20.0} # R_m = 10*10^6 Ohm
		
	d = 1.5 					# [ms] delay. same for all synapses
	
	# statistics constants	
	n_cells_per_vm = 2
	n_cells_to_raster = 30

	def build(self):
		self.inputs = Create("poisson_generator",self.n_minicolumns+1)
		self.input = Create(RC
		self.spikedetectors = Create("spike_detector", self.n_minicolumns+1)
		self.raster_spikedetector = Create("spike_detector")
		self.voltmeters = Create("voltmeter",self.n_minicolumns+1)
		
		#create populations of excitatory cells
		for n in range(self.n_minicolumns):
			mc = Create(HC.neuron_model, HC.N_E, {"t_ref":HC.t_ref})
			self.minicolumns.append(mc)
		
		#create population of inhibitory cells	
		self.basket_cells = Create(HC.neuron_model, HC.N_I*self.n_minicolumns, {"t_ref":5.0})
			
		# set rate for external poisson generators providing input to MC:s and BC:s 
		self.set_inrates(self.inrates)
		self.built = True
		
	def set_inrates(self, inrates):
		#set rate of external poisson input for respective minicolumn, last rate ev drive to basketcells
		for i,rate in enumerate(inrates):						#check way to avoid loop?
			SetStatus([self.inputs[i]],{"rate":rate})
			
	def connect(self):
		if not self.built:
			self.build()
		CopyModel("static_synapse", "excitatory", {"weight":HC.epsc_recurrent, "delay":HC.d})
		CopyModel("static_synapse", "inhibitory", {"weight":HC.ipsc_recurrent, "delay":HC.d})
		#CopyModel("static_synapse", "external", {"weight":HC.epsc_external, "delay":HC.d})
		
		#connect pyramidal cells internally and to basketcells, and basketcells to all pyramidal cells.
		#connect all ex populations to their spikedetector
		raster_neurons = []
		
		for i,mc in enumerate(self.minicolumns):
			RandomConvergentConnect(mc, mc, int(0.75*HC.N_E), model="excitatory")
			RandomConvergentConnect(mc, self.basket_cells, int(0.25*HC.N_E), model="excitatory")
			ConvergentConnect(self.basket_cells, mc, model="inhibitory")
			ConvergentConnect([self.inputs[i]], mc, model="excitatory")
			ConvergentConnect(mc[:40], [self.spikedetectors[i]])
			raster_neurons += mc[0:HC.n_cells_to_raster+1]
		
		raster_neurons +=self.basket_cells[0:HC.n_cells_to_raster+1]
		ConvergentConnect(raster_neurons, self.raster_spikedetector)
		
		#connect BC population to external input and spikedetector
		ConvergentConnect(self.inputs[-1:], self.basket_cells, model="excitatory")
		print(self.spikedetectors[-1:])
		ConvergentConnect(self.basket_cells, self.spikedetectors[-1:])
			
		#connect n cells from inh. and exc. populations to voltmeters
		for i,mc in enumerate(self.minicolumns):
			ex_cells = mc[0:HC.n_cells_per_vm]
			ConvergentConnect([self.voltmeters[i]], ex_cells)
		in_cells = self.basket_cells[0:HC.n_cells_per_vm]		
		ConvergentConnect([self.voltmeters[self.n_minicolumns]], in_cells)
			
		self.connected = True
			
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
		
		self.hc_average = sum(self.spike_statistics)/self.n_minicolumns
		
		#inhibitory spiking rate
		n_spikes = GetStatus([self.spikedetectors[-1]], "n_events")[0]
		spike_rate = 1000.0*n_spikes/(HC.t_sim*HC.N_I*self.n_minicolumns)
		self.spike_statistics.append(spike_rate)
		
		self.finished = True
	
	def simulation_info(self):	
		info = ""
		info+= "--------------------------------------\n"
		info+= "HYPERCOLUMN SIMULATION\t\t" + self.date_time + "\n"
		info+= "--------------------------------------\n"
		info+= "Number of minicolumns: \t\t%d\n" % self.n_minicolumns
		info+= "N_E per minicolumn:\t\t%d\n" %HC.N_E
		info+= "N_I per minicolumn:\t\t%d\n" % HC.N_I
		info+= "\n"
		info+= "Average spiking rate MC:s:\n"
		for i,spikerate in enumerate(self.spike_statistics[:-1]):
			info+= "\t\t\t\t%d:  %.2f Hz \t\t(%.2f)\n" % (i,spikerate,self.inrates[i])  
		info+="\n"
		info+="Average spiking rate HC:\t%d\n" % self.hc_average
		info+= "Average spiking inh. population\t\t%d\n" % self.spike_statistics[-1]
		info+= "\n"
		info+= "---------------------------------------\n\n"
		return info
		

	def voltage_traces(self, mcs=[]):
		"""Creates a figure and plots voltage_traces from given MC:s in each subplot"""
		if mcs == [] or len(mcs)>self.n_minicolumns:
			voltmeters = self.voltmeters
		else:
			voltmeters = []
			for i in mcs:
					voltmeters.append(self.voltmeters[i])
		subplot = len(voltmeters)*100+10
		plt.figure()
		for i, vm in enumerate(voltmeters[::-1]):
			subplot +=1	
			plt.subplot(subplot)
			simple_plots.voltage_trace([vm])
		#plt.show()	
		

	def hist(self,minicolumns):
		pass
	
	def __init__(self, n_minicolumns = 4, inrates=None):
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
		
		self.built = False
		self.connected = False
		ResetKernel()	
		
