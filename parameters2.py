params = {

	##############################################
	# PARAMETERS FOR HYPERCOLUMNS SIMULATION
	##############################################
	
	"t_sim":2000.0,				# [ms]
	"neuron_model":"iaf_neuron",
	
	# HYPERCOLUMN PARAMETERS
	"N_E":30,					# nr of pyramidal cells per MC
	"N_I_per_mc":4,				# nr of basket cells per MC, however modelled as an independent population size N_I*n_minicolumns
	"n_minicolumns":3,			# nr of minicolumns for HC

	"pyr_noise_rate":4000.0,    # [Hz] 4500 brings neurons to rheobase. 4000 weak sponateous firing
	"bas_noise_rate":0.0,		# [Hz] 5000 brings neuron to rheobase. 4500 weak spontaneous firing
	
	# SYNAPTIC WEIGHTS (CONTROL PSP size) 
	"epsc_pyr_pyr":120.0, 		# [pA] amplitude exc. postsynaptic current 
	"ipsc_bas_bas":-16.0,
	"epsc_pyr_bas": 16.0,			# [pA] amplitude exc from pyr to bas
	"epsc_ext_pyr":120.0,	
	"epsc_ext_bas":16.0,			#
	"epsc_noise_pyr":15.0,		# [pA] amplitude noise input. Gives EPSP = 0.1 mV
	"epsc_noise_bas":5.0,		# [pA] amplitude exc. postsynaptic current. 3.5 gives EPSP = 0.1 mV, 5.0 EPSP 0.14 mV
	"ipsc_bas_pyr":-340.0,		# [pA] amplitude inh. postsynaptic current
	
	"pyramidal_params":{"V_reset":0.0, "tau_m":20.0, "V_m":0.0, "E_L":0.0, \
				     "tau_syn":5.0, "t_ref":5.0, "C_m":1000.0, "V_th":20.0}, # R_m = 20*10^6 Ohm

	"basket_params":{"V_reset":0.0, "tau_m":10.0, "V_m":0.0, "E_L":0.0, \
				     "tau_syn":5.0, "t_ref":2.0, "C_m":200.0, "V_th":20.0},  # R_m = 50*10^6 Ohm
	# RANDOMIZED "SIZE" AND START MEMBRANE POTENTIAL
	"C_m_mean_bas":200.0,		
	"C_m_std_bas":1.0,			# biological value = 20?
	"C_m_mean_pyr":1000.0,	
	"C_m_std_pyr":1.0,			# biological value = 50?	
	"V_init_mean":5.0,
	"V_init_std":5.0,
	
	"synapse_std_fraction":0.25,
		
	"d":2.0,  					# [ms] delay. same for all synapses	
	
	# CONNECTION PROBABILITIES	
	"p_pyr_pyr":0.3,				# probability of ex.connection between pyr cells in same MC
	"p_bas_bas":0.0,
	"p_pyr_bas":0.7, 			# probability of ex.connection between pyr cells and basket cells
	"p_bas_pyr":0.7,	    		# probability of inh.connection between basket cells and pyr cells
	"p_in_pyr":0.5,				# probability of pyr cell getting exc. input
	"p_in_bas":0.4,				# input directly to bc. same epsc as from pyramidal
	
	# CONNECTION TYPES
	"mult":True,
	"aut":True,

	# EVALUATION PARAMETERS
	"n_cells_per_vm": 10,
	"n_cells_to_raster":30,

}
