import numpy as np
import pickle
import data_analysis as da
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import matplotlib.ticker as ticker
import plot_parameters as ppar

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def FIR_score_and_graph(all_inputs, outputs, n_th = 3, test="1234", z_max = 200, histogram = True  ):
	n_data, dim_data = len(all_inputs), len(outputs[0])
	n_mc = dim_data-1
	y_max = 250
	fcolor = {"1234":"#D6EBFF", "1234WTA":"m", "average":"0.92","1200":"#D6EBFF", "1200WTA":"m", "1333":"m", "1111":"1.0"}
	
	inhibitory = [x[-1] for x in outputs]
	outputs = [x[:-1] for x in outputs]
	averages_in = [np.sum(x)/n_mc for x in all_inputs]
	averages = [np.sum(a)/n_mc for a in outputs]

	norm_activities = np.zeros((n_mc,n_data))
	for i in range(n_mc):
		norm_activities[i] = [x[i]/(averages[j]*n_mc) for j,x in enumerate(outputs)]

	# bar diagrams
	outputs_to_bars = [x for i,x in enumerate(outputs) if i%n_th==0]
	bars = np.array(outputs_to_bars).flatten()
	# colors
	bar_colors = ['#00FFFF','#3333CC','#CC0099']
	
	if len(bar_colors)> n_mc:
		bar_colors.remove()
	extra_colors = np.linspace(0.75,0.8,n_mc-len(bar_colors))
	for extra in extra_colors:
		bar_colors.append(str(extra))

	bar_colors_all = bar_colors*(len(outputs_to_bars)/n_th)

	x = np.arange(0,len(bars))
	dots = [y+n_mc/4 for j,y in enumerate(x) if j%n_mc==0]
	averages_plot = [a for i,a in enumerate(averages) if i%n_th==0]
	inhibitory_plot = [d[-1] for i,d in enumerate(outputs) if i%n_th==0]

	best_seq, start, stop, best_q = longest_stable_average(np.array([sum(g)/len(g) for g in all_inputs]),\
															outputs,test, False)
	print "Best sequence:", best_seq, "start:",start, "stop:", stop, "best_q", best_q

	if histogram:
		fig, ax = plt.subplots()
		for j in range(len(norm_activities)):
			norm_activity_plot = [n*100 for i,n in enumerate(norm_activities[j]) if i%n_th==0]#
			plt.plot(dots,norm_activity_plot,marker='o', mew=0.1,linewidth=0.75, color=bar_colors[j])
		inhibitory_plot = [inh for i,inh in enumerate(inhibitory) if i%n_th==0]
		if best_q == 1:
			start, stop = 0, 0.5
		else:
			start, stop = best_seq[0], best_seq[1]
		plt.axvspan(start*n_mc/float(n_th), stop*n_mc/float(n_th)+n_mc, facecolor=fcolor[test], linewidth=0.5)
		plt.bar(x,bars,color=bar_colors_all, linewidth=0.1, width=1.0)
		plt.plot(dots,averages_plot,linewidth=1.5, color ="k")#"#4D4D4D"
		l, = plt.plot(dots,inhibitory_plot,"--",color="0.1", linewidth=0.75)
		l.set_dashes([3,4])
		print np.round(dots,-2)

		plt.ylim([0,y_max])
		ax.set_xticks(np.round(np.linspace(0,len(bars)-3,5),0))
		ax.set_xticklabels([50,250,500,2500,10000])
		ax.set_yticks([0,50,100,150,200])
		plt.tick_params(axis="y", which="both", right = "on", left = "on", labelleft="on", labelright="off")
		plt.tick_params(axis="x", which="both", top="off", bottom="on")

		if n_th == 2:
			plt.xlim([0,len(all_inputs)*n_mc/n_th+2])
		else:
			plt.xlim([0,len(all_inputs)*n_mc/n_th])

		plt.xlabel("Average input (Hz)")
		plt.ylabel("Output (Hz)")
		ax.xaxis.labelpad = 2
		ax.yaxis.labelpad = 2

		ax2 = ax.twinx()
		plt.ylim([0,y_max])
		#ax2.set_yticklabels([0,20,40,60,80,100])
		#ax2.yaxis.set_ticks([0,20,40,60,80,100])
		ax2.set_yticklabels([0,50,100])
		ax2.yaxis.set_ticks([0,50,100])
		plt.tick_params(axis="y", which="both", right = "on", left = "off", labelleft="off", labelright="on")
		#	ax2.set_yticklabels([])
	#	ax.set_yticklabels([])

		text = "FIR score:  " + "%0.2f" %best_q
		ax.text(0.05, 0.8, text, transform = ax.transAxes, fontsize=12)

		#plt.savefig("fir_test.svg")
		plt.show()

	print "Return utilility"
	return best_seq, start, stop, best_q

def surface_3d(X1,X2,Z, z_lim = 100, x1_label = "Input mc1 (Hz)", x2_label="Input mc2 (Hz)", z_label="Average output (Hz)"):
    X1,X2 = np.meshgrid(X1,X2)

    fig = plt.figure()
    ppar.set2perpage3D()
    ax = fig.gca(projection="3d")
    ax.set_zlim(0,100)
    ppar.ax_stuff2(ax)

    surf = ax.plot_surface(X1, X2, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0.0, \
    antialiased=False) 

    ax.zaxis.set_major_locator(LinearLocator(6))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.f'))
    xLabel = ax.set_xlabel(x1_label)
    yLabel = ax.set_ylabel(x2_label)
    zLabel = ax.set_zlabel(z_label)
    ax.xaxis._axinfo['label']['space_factor'] = 1.9
    ax.yaxis._axinfo['label']['space_factor'] = 2.0
    ax.zaxis._axinfo['label']['space_factor'] = 1.7
    ax.dist = 9

    cb = fig.colorbar(surf, shrink=0.5, aspect=15, pad = 0.05)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()

    plt.show()

def contour_2d(X1,X2, Z):
    ppar.set2perpage2D()
    fig, ax = plt.subplots(1)
    ppar.ax_stuff(ax)
    clevels = [0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.85, 5.0, 10.0, 20.0, 50.0, 100.0]

    fmt = '%0.1f'
    cs = plt.contour(X1,X2,Z,clevels, norm=LogNorm(), cmap=cm.ocean)
    plt.clabel(cs, [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0], inline=True, fmt=fmt, fontsize=8,textcolor="k", norm = LogNorm())

   # plt.colorbar(cs)
    plt.ylabel("Input mc 1 (Hz)")
    plt.xlabel("Input mc 2 (Hz)")
    plt.show()



###############################################################################################
# Function below calculates FIR test score from outputs
###############################################################################################

test_params = {"1333":{"p":[10.0, 10.0, 10.0, 10.0, 0.5],"relations":[1.25, 0.0, 0.0]}, 
			  "1234WTA":{"p":[1000.0, 1000.0, 1000.0, 1000.0, 0.2],"relations":[0.0, 1.1, 0.0]},
			   "1234":{"p":[10.0, 0.25, 0.25, 0.25, 0.20], "relations":[1.50,1.25,1.15]},
			   "1200":{"p":[0.5,0.25,1000.0,1000.0,0.20],"relations":[1.5, 0.0, 0.0],"min_relations":[0.0, 0.0] },
			   "1200WTA":{"p":[1000.0,1000.0,1000.0,1000.0,0.20],"relations":[4.0, 0.0, 0.0]},
			   "1111":{"p":[0.25,0.25,0.25,0.25,0.20],"relations":[0.0, 0.0, 0.0]},
		       "average":{"p":[1000.0,1000.0,1000.0,1000.0, 0.20],"relations":[0.01, 0.01, 0.01] }}

def longest_stable_average(input_averages, statistics, test="1234", debug = False):
	
	n_data = len(input_averages)
	n_mc = len(statistics[0])
	ps = test_params[test]["p"]
	relations = test_params[test]["relations"]
	
	averages = [sum(m)/float(n_mc) for m in statistics]	
	outputs = np.zeros((n_mc,n_data))
	np.set_printoptions(precision=3)
# 
	# change matrix dimensions 
	for i in range(n_mc):
		outputs[i] = [x[i:i+1]/(averages[j]*n_mc) if x[i:i+1]/(averages[j]*n_mc)>0.01 else 0.01 for j,x in enumerate(statistics)]
		
	s_averages = np.zeros((n_data, n_data))
	relation_averages = np.zeros((n_mc, n_data, n_data))
	
	# calculate average for each possible sequence O(n^2)
	for i in range(n_data):
		for j in range(i,n_data):
			s_averages[i][j] = sum(averages[i:j+1])/len(averages[i:j+1])
			
			for k in range(n_mc):
				relation_averages[k][i][j] = sum(outputs[k][i:j+1])/len(outputs[k][i:j+1])
	
	# calculate if the sequences are withing +/- p perecent of their average 
	s_averages_ok = np.zeros((n_data, n_data))
	mc_s_averages_ok = np.zeros((n_mc, n_data, n_data))
	
	for i in range(n_data):
		for j in range(i,n_data):
			if allowed_sequence(averages[i:j+1],s_averages[i][j],ps[-1]):
				s_averages_ok[i][j] = 1
			for k in range(n_mc):
				if allowed_sequence(outputs[k][i:j+1], relation_averages[k][i][j], ps[k] ): 
					mc_s_averages_ok[k][i][j]=1

	# calculate if the relations between mc:s are ok for sequences
	mc_relations = np.zeros((n_mc-1,n_data,n_data))

	for k in range(n_mc-1):
		for i in range(n_data):
			for j in range(i,n_data):
				if mc_relations_ok(outputs[k][i:j+1], outputs[k+1][i:j+1], relations[k]):
					#print "mc_relations", mc_relations
					mc_relations[k][i][j]=1
					
				if "min_relations" in test_params[test].keys() and mc_relations[k][i][j] == 1:
					mc_relations[k][i][j] = absolute_relations_ok( outputs[k][i:j+1], test_params[test]["min_relations"])
					
					
	# get best allowed sequence, based on max quota between the average input for min and max
	best_seq = (-1,-1)
	best_q = 1
	for i in range(n_data):
		for j in range(i,n_data):
			if s_averages_ok[i][j]:
				ok = True
				for k in range(n_mc-1):						
					if not (mc_s_averages_ok[k][i][j] and mc_relations[k][i][j] ):
						ok = False
				if not mc_s_averages_ok[n_mc-1][i][j]:
						ok = False
				if ok:				
					input_q = input_averages[j]/input_averages[i]
					#if (j-i)>best_q:
					if input_q > best_q:
						best_q = input_q
						best_seq = (i,j)
	
	if debug:	
		print "Statistics:", statistics
		print "n_mc",n_mc
		#print "averages"
		#print averages	
		print "s_averages"
		print s_averages
		print "relation averages"
		for relations in relation_averages[0:2]:
			print relations			
		print "s_averages_ok"
		print s_averages_ok
		print "relation averages_ok"
		for rel_av in mc_s_averages_ok[0:2]:
			print rel_av
		print "mc_relations_ok"
		for mc_r in mc_relations:
			print "----------------------------------------------------------------------------------"
			print mc_r
	
	return best_seq, input_averages[best_seq[0]], input_averages[best_seq[1]], best_q

def mc_relations_ok(outputs1, outputs2,pp):
	#print outputs1
	#print outputs2
	for i,o in enumerate(outputs1):
		if outputs2[i]/outputs1[i] < pp:
			return False
	return True
	
def absolute_relations_ok(vector,lim):
	return True
	#print "vector", vector
	#for o in vector:
		#print "o", o
		#if o < lim:
			#return False
	#return True
	
def allowed_sequence(data, average, pp):
	#print "average", average
	max_ok = (1.0+pp)*average
	min_ok = (1.0-pp)*average
	#print "data", data
	for d in data:
	#	print "d", d
		if d < min_ok or d>max_ok:
			return False
	return True
	
def divide_without_infinity(Y1,Y2):
	res = np.zeros((len(Y1),len(Y2[0])))
	for row in range(len(res)):
		for col in range(len(res[0])):
			if Y1[row][col] == 0 or Y2[row][col]/Y1[row][col]>20.0:
				res[row][col]=22
			elif Y2[row][col]/Y1[row][col]<0.05:
				res[row][col]=0.049

			else:
				res[row][col]=Y2[row][col]/Y1[row][col]
	return res


def save_obj(data, file_name):
	myfile = open( file_name + ".pkl", 'w')
	pickle.dump(data, myfile)

def load_obj(file_name ):
    with open(file_name + '.pkl', 'r') as f:
        return pickle.load(f)
        
 
