import nest
from nest import simple_plots
from nest import raster_plot
from matplotlib import pyplot as plt
import numpy as np
import hc_utils as hcu


class Visualization:

	def __init__(self, hc):
		self.hc = hc
		
	def visualize(self, plots = "hvrs"):
		"""2*2 grid with of 5*1 rows"""
		n_mcs = self.hc.params["n_minicolumns"]
		t_sim = self.hc.params["t_sim"]
		n_bins = 100	
		bins = np.linspace(0,t_sim,n_bins+1)
		bin_size =t_sim/n_bins
	
		# Spike histograms
		if "h" in plots:
			f, axarr= plt.subplots(n_mcs+1,sharex=True)
			plt.xlabel("Time (ms)")
			f.subplots_adjust(hspace=0.05)
			for i, sd in enumerate(self.hc.spikedetectors[::-1]):
				axarr[i].set_yticks([0,100])
				axarr[i].set_ylim([0, 200])
				axarr[i].set_xlim([0, t_sim])
				axarr[i].xaxis.labelpad = 2
				axarr[i].yaxis.labelpad = 2
				
				events = nest.GetStatus([sd],"events")[0]
				times = events["times"]
	
				#self.setTickSize(axarr[i])

		
				if not i == 0: 							# Om inte inhibitory pop..
					norm_factor = 1000.0/(bin_size*self.hc.params["N_E"]) 
				else: 											# ..som har annat antal neuroner
					norm_factor = 1000.0/(bin_size*self.hc.params["N_I"]) 
				h_weights = [norm_factor]*len(times)
				
				text2 = "%.1f" %(self.hc.spike_statistics[n_mcs-i])
				#axarr[i].text(0.85, 0.78, text2, transform = axarr[i].transAxes, fontsize = 10)
			
				if i == 0:
					axarr[i].set_ylim([0, 250])
					axarr[i].set_yticks([0,100,200])
					text =  "(" + str(int(sum(self.hc.inrates)*self.hc.params["p_in_bas"])) +")"
					#axarr[i].axhspan(0,t_sim, facecolor="grey", alpha =0.3)
					axarr[i].set_axis_bgcolor("0.90")
					#axarr[i].xaxis.set_ticks_position("bottom")
					#axarr[i].yaxis.set_ticks_position("left")
				else:
					text = "(" + str(int(self.hc.inrates[n_mcs-i])) + ")"
				axarr[i].yaxis.set_ticks_position("left")
				axarr[i].text(1.02, 0.20, text + "\n" + text2, transform = axarr[i].transAxes, fontsize=8)
				if len(times):
					axarr[i].hist(times, bins, weights=h_weights, linewidth=0.5)		# vikta staplarna for att fa i <Hz>			
			f.text(0.06, 0.5, 'Spike frequency (Hz)', ha='center', va='center', rotation='vertical', fontsize=8)
#		
		# Voltage trace plots	
		if "v" in plots:	
			f, axarr= plt.subplots(n_mcs+1,sharex=True)
			plt.xlabel("Time (ms)")
			f.subplots_adjust(hspace=0.05)
			
			for i, vm in enumerate(self.hc.voltmeters[::-1]):
				if i == 0:
					axarr[i].set_axis_bgcolor("0.90")
					text =  "(" + str(int(sum(self.hc.inrates)*self.hc.params["p_in_bas"])) +")"
				else:
					text = "(" + str(int(self.hc.inrates[n_mcs-i])) + ")"
				axarr[i].set_ylim([-5, 22])
				axarr[i].set_xlim([0, t_sim])
				axarr[i].set_yticks([0,10])
				axarr[i].xaxis.labelpad = 2
				axarr[i].yaxis.labelpad = 2
				#self.setTickSize(axarr[i],10)
				simple_plots.voltage_trace2(axarr[i], [vm])
				axarr[i].yaxis.set_ticks_position("left")
			
				axarr[i].text(1.02, 0.4, text, transform = axarr[i].transAxes, fontsize=8)	
				
			
			f.text(0.06, 0.5, 'Membrane voltage (mV)', ha='center', va='center', rotation='vertical', fontsize=8)
	
		#	plt.show()

		# Raster plots
		if "r" in plots:	
			f, ax = plt.subplots(1)
			ax.xaxis.labelpad = 1
			plt.xlabel("Time (ms)")
			
			all_text = ""
			for i in range(self.hc.params["n_minicolumns"]+1):
				if i==0:
					text =  "(" + str(int(sum(self.hc.inrates)*self.hc.params["p_in_bas"])) +")"
					ax.text(1.02, 0.92, text , transform = ax.transAxes, fontsize=8)
				else:
					text = "(" + str(int(self.hc.inrates[n_mcs-i])) + ")"+"\n"
					all_text += text +"\n\n" 
			
			ax.text(1.02, -0.1, all_text , transform = ax.transAxes, fontsize=8)
	
			ax.set_yticks(np.arange(self.hc.all_pyr[0], self.hc.all_pyr[-1]+1, 30))
			ax.grid(True,'major','y',linestyle='-')
			plt.setp(ax.get_yticklabels() , visible=False)
	
			raster_plot.from_device(self.hc.raster_spikedetector, title="")
			plt.ylim(self.hc.all_pyr[0],self.hc.basket_cells[-1]+0.5)
			#ax.axhspan(self.hc.basket_cells[-1],self.hc.basket_cells[0]-0.5, facecolor="grey", alpha =0.3)
			ax.axhspan(self.hc.all_pyr[-1]+0.5,self.hc.basket_cells[-1]+0.5, facecolor="0.9", linewidth=0.5)
			
			plt.ylabel("Model neuron id")
		
		# Plot over spike distribution
		if "s" in plots:
			
			f, axarr= plt.subplots(n_mcs+1,sharex=False)
			f.subplots_adjust(hspace=0.1)
			
			for m in range(n_mcs+1):
			
				spike_statistics_mc, neurons = hcu.spike_frequency_analysis(self.hc, n_mcs-m)
				std, mean = np.std(spike_statistics_mc), np.mean(spike_statistics_mc)
				spike_statistics_mc.sort()
		
				axarr[m].bar(range(len(spike_statistics_mc)),spike_statistics_mc,linewidth=0.5, width=1.0)
				l, = axarr[m].plot([0,30],[mean, mean],"--k",linewidth=0.75)
				l.set_dashes([3,4])
				
				#text = "Spikefrequency: %.2f +/- %.2f" % (mean,std)
				#axarr[m].text(0.5, 0.7, text, transform = axarr[m].transAxes, fontsize=8)
				plt.setp(axarr[m].get_xticklabels() , visible=False)
				axarr[m].yaxis.set_ticks_position("left")
				axarr[m].xaxis.set_ticks_position("none")
				axarr[m].set_xlim([0,30])
				if m != 0:
					axarr[m].set_ylim([0, 200])
					axarr[m].set_yticks([0,100])
					text = "(" + str(int(self.hc.inrates[n_mcs-m])) + ")"
				else:
					axarr[m].set_axis_bgcolor("0.90")
					axarr[m].set_ylim([0, 200])
					axarr[m].set_yticks([0,100])
					text =  "(" + str(int(sum(self.hc.inrates)*self.hc.params["p_in_bas"])) +")"
				axarr[m].text(1.02, 0.4, text, transform = axarr[m].transAxes, fontsize=8)	
				plt.xlabel("Model neuron id")
		plt.show()		


		
	def setTickSize(self, ax, fs=10):
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
			item.set_fontsize(fs)
		#for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		#	label.set_fontsize(8)
			
		
		
