params = {

	##############################################
	# PARAMETERS FOR HYPERCOLUMNS SIMULATION
	# iaf_neuron with exponential psc, "iaf_psc_exp"
	# parameters "normal" not tweaked to give gain change
	# that is represents Model A
	##############################################
	
	"info":"Parameters for Model A1",
	"t_sim":500.0,				# [ms]
	
	# NEURON PARAMETERS	
	"neuron_model":"iaf_cond_exp",
	"pyramidal_params":{"V_reset":0.0, "g_L":5.2 ,"V_m":0.0, "E_L":0.0, "E_ex":120.0, "E_in":-10.0, \
				     "tau_syn_ex":6.0, "tau_syn_in":6.0, "t_ref":3.5, "C_m":70.0, "V_th":15.0},  
					#-> tau_mem = 13.5 ms, R_m = 190 MOhm
 
	"basket_params":{"V_reset":0.0, "g_L":0.56, "V_m":0.0, "E_L":0.0, "E_ex":120.0, "E_in":-10.0, \
				     "tau_syn_ex":6.0, "tau_syn_in":6.0, "t_ref":2.0, "C_m":7.5, "V_th":15.0}, 	
					#-> tau_mem = 13.5 ms, R_m = 1800 MOhm
	
	# HYPERCOLUMN PARAMETERS
	"N_E":30,					# nr of pyramidal cells per MC
	"N_I":16,				# nr of basket cells per MC, however modelled as an independent population size N_I*n_minicolumns
	"n_minicolumns":4,			# nr of minicolumns for HC

	"pyr_noise_rate":5200.0,    	# [Hz] 5200 brings neurons appr. to rheobase. 4000 weak sponateous firing
	"bas_noise_rate":5000.0, # [Hz] 5500 brings neuron to rheobase. 4500 weak spontaneous firing, 		1250
	
	# SYNAPTIC WEIGHTS (CONTROL PSP size)
	"epsp_ext_pyr":0.17,		# 0.17 peak conductance in [nS] 
	"epsp_pyr_pyr":0.17, 		# 0.17 
 	"epsp_noise_pyr":0.020,		# 0.020
	"ipsp_bas_pyr": - 40.0,	# -2.6:-1.2mV     -12.0

	"epsp_ext_bas": 0.001,	# 0.009		0.03
	"epsp_pyr_bas": 0.001,	# 0.009		0.004
	"epsp_noise_bas":0.0020,	# 0.0020 5.35 + 4.0   
	"ipsp_bas_bas": - 0.12,	#NO
	
	# RANDOMIZED "SIZE" AND START MEMBRANE POTENTIAL used only for randomized versions
	"C_m_mean_pyr": 70.0,		
	"C_m_sd_fraction_pyr": 0.1,			# biological value?	
	"C_m_mean_bas":7.5,		
	"C_m_sd_fraction_bas": 0.1,			# biological value?
	"V_init_mean":  5.0,
	"V_init_sd_fraction":   1.0,
	#"synapse_std_fraction":0.1,
	"Cm_Vm_normal_clip":1.0,
	
	"d":1.5,  					# [ms] delay. same for all synapses	
	
	# CONNECTION PROBABILITIES	
	"p_pyr_pyr":0.2, 			# probability of ex.connection between pyr cells in same MC
	"p_bas_bas":0.0,
	"p_pyr_bas":0.7,		# probability of ex.connection between pyr cells and basket cells
	"p_bas_pyr":0.7,    		# probability of inh.connection between basket cells and pyr cells
	"p_mc_mc":0.0,				# (lateral) connection probability between mc:s
	"p_in_pyr":None,			# probability of pyr cell getting exc. input
	"p_in_bas":0.05,			    # input directly to bc. same epsp as from pyramidal
	
	"input_std_fraction": 0.10,			# for randominzing the amount of input to mc
	"input_std_fraction_bas":0.10,		# and basket cells	
	"input_normal_clip":2.0,
	
	# CONNECTION TYPES
	"mult":True,
	"aut":True,
	"direct_inhibition":False,

	# EVALUATION PARAMETERS
	"n_cells_per_vm": 5,
	"n_cells_to_raster":30,

	"random_seeds":[2,7,10,15,16]   # 2 --> 5.35, 7 --> 4.05, 10 --> 4.65
}
