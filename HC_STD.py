from nest import *
import numpy as np
from scipy.stats import truncnorm
from HC_ModelA import HC

import random
import sys
import time

#######################################################################
# Inherits from HC_ModelA but defines new build and connect methods 
# for including short term depression and a population of "input neurons"
#######################################################################


class HC_STD(HC):
	"""Hypercolumn class with variable number of minicolumns"""

	rates_other = [100.0, 400.0, 700.0, 1000.0, 1300.0] 
	rates = np.concatenate((np.arange(0, 3001, 500.0, "float"),np.arange(2000, 2000, 250, "float")),0)
	
	
	def build(self):		
		# Each neuron will have its own input rate, and therefore an indivdual Poisson generator
		self.inputs = Create("poisson_generator",self.params["n_minicolumns"]+1)
		if self.params["direct_inhibition"]:
			self.direct_inhibition = Create("poisson_generator")
		self.inputs[-1] = Create("poisson_generator",self.params["N_I"])
		self.pyr_noise_input = Create("poisson_generator")				#noise input to all exc neurons
		self.bas_noise_input = Create("poisson_generator")				#noise input to all inh neurons
		self.spikedetectors = Create("spike_detector", self.params["n_minicolumns"]+1)
		self.l4_spikedetectors = Create("spike_detector", self.params["n_minicolumns"])
	
		self.raster_spikedetector = Create("spike_detector")
		self.voltmeters = Create("voltmeter",self.params["n_minicolumns"]+1)
		self.all_pyr = []												#list of all gids for pyramidal cells, for convenience 
		
		#layer4 
		for n in range(self.params["n_minicolumns"]+1):
			l4 = Create("parrot_neuron", self.params["N_E_4"])
			self.layer4.append(l4)
		
		#create populations of excitatory cells
		for n in range(self.params["n_minicolumns"]):
			mc = Create(self.params["neuron_model"], self.params["N_E"],\
				        params=self.params["pyramidal_params"])
			self.all_pyr = self.all_pyr + mc  	
			self.minicolumns.append(mc)
	
		#create population of inhibitory cells	
		self.basket_cells = Create(self.params["neuron_model"], self.params["N_I"], \
								   params=self.params["basket_params"])		
		self.randomize_network_and_inputs() 	
	
		self.built = True	

	def connect(self):
		if not self.built:
			self.build()
		
		# depressing synapse model
		tsodyks_params = self.params["tsodyks_params"]
		tsodyks_params["weight"] = self.params["epsp_ext_pyr"]*2
		tsodyks_params["delay"] = self.params["d"]
		syn = CopyModel("tsodyks2_synapse", "pyr_std", params = tsodyks_params)
	
		CopyModel("static_synapse", "pyr_recurrent", {"weight":self.params["epsp_pyr_pyr"], "delay":self.params["d"]}) 
		CopyModel("static_synapse", "bas_external", {"weight":self.params["epsp_ext_bas"],"delay":self.params["d"]})
		CopyModel("static_synapse", "bas_recurrent", {"weight":self.params["ipsp_bas_bas"],"delay":self.params["d"]})	
		CopyModel("static_synapse", "pyr_inhibitory", {"weight":self.params["ipsp_bas_pyr"],"delay":self.params["d"]})
		CopyModel("static_synapse", "bas_excitatory", {"weight":self.params["epsp_pyr_bas"],"delay":self.params["d"]})
	
		# Choise between fixed number of incoming or outgoing connections
		if self.params["mult"] == True:
			connect = RandomDivergentConnect
		else:
			connect = RandomConvergentConnect
	
		# Connections..	
		self.raster_neurons = []
		for i,mc in enumerate(self.minicolumns):	
			DivergentConnect([self.inputs[i]],self.layer4[i])
			connect(self.layer4[i],  mc, int(self.params["p_in_pyr"]*self.params["N_E_4"]), \
					model ="pyr_std", options={'allow_multapses':self.params["mult"], 'allow_autapses':self.params["aut"]})
			
			connect(self.layer4[i],  self.basket_cells, int(self.params["p_in_bas"]*self.params["p_in_pyr"]*\
					self.params["N_E_4"]), model ="bas_external", options={'allow_multapses':self.params["mult"], 'allow_autapses':self.params["aut"]})
			
			connect(mc, mc, int(self.params["p_pyr_pyr"]*self.params["N_E"]), \
					model="pyr_recurrent", options={'allow_multapses':self.params["mult"], 'allow_autapses':self.params["aut"]})
			connect(mc, self.basket_cells, int(self.params["p_pyr_bas"]*self.params["N_E"]), \
					model="bas_excitatory", options={'allow_multapses':self.params["mult"]})
			connect(self.basket_cells, mc, int(self.params["p_bas_pyr"]*self.params["N_I"]), \
					model="pyr_inhibitory", options={'allow_multapses':self.params["mult"]})			
			DivergentConnect(self.pyr_noise_input, mc, self.params["epsp_noise_pyr"], self.params["d"])
	

			# .. to device
			self.raster_neurons += mc[0:self.params["n_cells_to_raster"]+1]
	
		if self.params["p_mc_mc"]>0:
			self.mc_mc_connect()
		
		# .. external noise bc pop, recurrent bc pop
		DivergentConnect(self.bas_noise_input, self.basket_cells, self.params["epsp_noise_bas"], self.params["d"])
		RandomDivergentConnect(self.basket_cells, self.basket_cells, int(self.params["p_bas_bas"]*self.params["N_I"]),\
					model="bas_recurrent")
		self.raster_neurons +=self.basket_cells[0:self.params["n_cells_to_raster"]+1]
		
		self.connect_devices()
		self.connected = True		
	
	
	def set_inrates(self, inrates, inh_rate = None):	
		""" Set rates for excitatory and inhibitory input poisson generators. Inhibitory
			rates calculated from excitatory rates if inh_rate equals None"""
		inrates = np.array(inrates)/(self.params["N_E_4"]*self.params["p_in_pyr"])
		
		if len(inrates) != self.params["n_minicolumns"]:
			raise RuntimeError("n_inputs does not equal n_minicolumns")
		self.inrates = inrates
		print "Inrates:", inrates
		
		# Input to pyramidal cells
		for i,rate in enumerate(self.inrates):
			SetStatus([self.inputs[i]], "rate", float(rate))
		
		SetStatus(self.pyr_noise_input,{"rate":self.params["pyr_noise_rate"]})
		SetStatus(self.bas_noise_input,{"rate":self.params["bas_noise_rate"]})
		
	 
	def __init__(self, params):		
		ResetKernel()			
		self.params = params				# parameter dictionary, test script reads this from file
		self.minicolumns = []				# list of nodeslists. MCs added here during build	
		self.layer4 = []
		self.inrates = None			
		self.setRandomSeeds(params["random_seeds"][1])	
		SetKernelStatus({"resolution":0.1})
		self.finished, self.built, self.connected =  False, False, False
