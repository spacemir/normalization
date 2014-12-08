import nest
from matplotlib import pyplot
import numpy
import array

def activity_trace2(device, sim_time, binsize=10):
	"""binsize and simtime in ms"""
	n_events = nest.GetStatus(device, "n_events")[0]
	ev = nest.GetStatus(device, "events")[0]
	event_times = ev["times"]
	
	n_bins = int(sim_time/binsize)
	hist_bins = numpy.zeros(n_bins)
		
	for t in event_times:
		h_bin = int(t//10)
		hist_bins[h_bin]+=1
	
	hz_trace = hist_bins/(binsize*1e-3)		
	
	return hz_trace, hist_bins 	

def from_memory(device):
	""" returns two dictionaries with times and voltages for all neurons"""
	
	ev = nest.GetStatus(device, "events")[0]
	potential = ev["V_m"]
	senders = ev["senders"]
	v= {}
	t= {}
	times = ev["times"]
	for s, currentsender in enumerate(senders):
		if not v.has_key(currentsender):
			v[currentsender] = array.array('f')
			t[currentsender] = array.array('f')
		v[currentsender].append(float(potential[s]))
		t[currentsender].append(float(times[s]))
	return t, v	
		

def voltage_trace(ax, device, standalone = False, n_neurons = None):
	times, voltages = from_memory(device)
	keys = voltages.keys()
	neurons = keys
	if n_neurons == None or n_neurons>len(keys):
		neurons=keys
	else:
		n_neurons = keys[0:n_neurons+1]
	time_values = numpy.array(times[neurons[0]])
	for n in neurons:
		ax.plot(time_values, voltages[n])
	#pyplot.title("Membrane potential [mV]") 

	return None
	

def voltage_trace2(ax, device, n_neurons = None):
	times,voltages = from_memory(device)
	keys = voltages.keys()
	neurons = keys
	
	if n_neurons == None or n_neurons>len(keys):
		neurons=keys
	else:
		n_neurons = keys[0:n_neurons+1]
	time_values = numpy.array(times[neurons[0]])
	for n in neurons:
		ax.plot(time_values, voltages[n])
	#pyplot.title("Membrane potential [mV]") 
	return ax
	

def histogram(device, hist_binwidth=5.0, standalone = False):
	ev = nest.GetStatus(device, "events")[0]
	

	
	
