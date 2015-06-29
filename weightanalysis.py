import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
import numpy as np
import zen
import random
import os
import time
import seaborn as sns

def load_network(networkname):
	G = zen.io.edgelist.read(networkname, directed=True, weighted=True)
	return G

def degree_distribution(G):
	in_degrees = [G.in_degree_(nidx) for nidx in G.nodes_iter_()]
	out_degrees = [G.out_degree_(nidx) for nidx in G.nodes_iter_()]
	degrees = [G.degree_(nidx) for nidx in G.nodes_iter_()]
	return (in_degrees, out_degrees, degrees)

def plot_deg_distrib(G):
	(in_deg, out_deg, deg) = degree_distribution(G)
	plt.figure(1)
	plt.subplot(131)
	plt.hist(in_deg)
	plt.subplot(132)
	plt.hist(out_deg)
	plt.subplot(133)
	plt.hist(deg)
	plt.show()

def get_controlled_nodes(G):
	matched_edges = zen.maximum_matching_(G)
	nodes = G.nodes_()
	for eidx in matched_edges:
		if G.tgt_(eidx) in nodes:
			nodes = np.delete(nodes, np.where(nodes==G.tgt_(eidx)))
	return nodes, matched_edges

TYPE_SOURCE = 'src'
TYPE_EXTERNAL_DILATION = 'ext'
TYPE_INTERNAL_DILATION = 'int'

def control_profile_(G, controls, matching):
	types = dict()
	num_source = 0
	num_external_dilations = 0
	num_internal_dilations = 0
	sources = [G.src_(e) for e in matching]
	targets = [G.tgt_(e) for e in matching]
	for i, control in enumerate(controls):
		end = 0
		if G.in_degree_(control) == 0: # control is a source
			types[control] = TYPE_SOURCE
			num_source += 1
		else:
			if control in sources:
				idx = sources.index(control)
				while targets[idx] in sources: # find the end of the stem
					idx = sources.index(targets[idx])
				end = targets[idx]
			else:
				end = control

			if G.out_degree_(end) == 0: # end of a stem is a sink
				types[control] = TYPE_EXTERNAL_DILATION
				num_external_dilations += 1
			else:
				types[control] = TYPE_EXTERNAL_DILATION
				num_internal_dilations += 1

	return (float(num_source)/len(controls), float(num_external_dilations)/len(controls), float(num_internal_dilations)/len(controls)), types

def find_one_min_weight(G, w):
	min_weights = []
	abs_min_weight = min(w)
	for eidx, wei in G.edges_iter_(weight=True):
		if wei >= abs_min_weight-10**(-6) and wei <= abs_min_weight+10**(-6):
			min_weights.append([eidx, wei])
	chosen_one = min_weights[random.randint(0, len(min_weights)-1)]
	return chosen_one

def find_one_max_weight(G, w):
	max_weights = []
	abs_max_weight = max(w)
	for eidx, wei in G.edges_iter_(weight=True):
		if wei >= abs_max_weight-10**(-6) and wei <= abs_max_weight+10**(-6):
			max_weights.append([eidx, wei])
	chosen_one = max_weights[random.randint(0, len(max_weights)-1)]
	return chosen_one

def avg_reducer(H, step, runs=100, rand=False, asc=True, filename=""):
	iter_number = H.size()/step + 1
	avg_controls = np.zeros((iter_number, 4))
	for i in range(0, runs):
		print 'Run %d' % i
		[controls, diff_controls, weight_cuts, node_removed, controlled_nodes] = reducer(H, step, rand, asc)
		avg_controls = avg_controls+controls/100
	if filename != "":
		np.savetxt("controls/"+filename+".txt", avg_controls)
	return avg_controls

def reducer(H, step, rand=False, asc=True):
	total_running_time = time.time()
	G = H.copy()
	w = [w for u,v,w in G.edges_iter(weight=True)]
	iter_number = G.size()/step + 1
	controls = np.zeros((iter_number, 4))	
	diff_controls = np.zeros(iter_number-1)
	weight_cuts = np.zeros(iter_number-1)
	controlled_nodes = np.ones((len(G), iter_number))
	remove_disconnected_nodes = False
	node_removed = np.zeros(iter_number-1)
	initial_size = len(G)

	for i in range(0, iter_number):
		print "Slice %d out of %d" % (i, iter_number)
		if i > 0: # First slice is the whole graph
			start_time = time.time()
			for j in range(0, step): # slice 'step' edges from the graph
				if G.size() > 0:
					chosen_one = ()
					if not rand:
						if asc:
							chosen_one = find_one_min_weight(G,w)
						else:
							chosen_one = find_one_max_weight(G,w)
					else:
						rem = random.randint(0,G.size()-1)
						eidx = G.edges_()[rem]
						wei = G.weight_(eidx)
						chosen_one = [eidx, wei]
					eidx = chosen_one[0]
					src = G.src_(eidx)
					tgt = G.tgt_(eidx)
					index = np.where(G.edges_()==eidx)[0][0]
					G.rm_edge_(eidx)
					w.pop(index)
					weight_cuts[i-1] = chosen_one[1]
					if remove_disconnected_nodes:
						if G.degree_(src) == 0:
							G.rm_node_(src)
							diff_node_removed[i-1] += 1
						if G.degree_(tgt) == 0:
							G.rm_node_(tgt)
							diff_node_removed[i-1] += 1
			node_removed[i-1] = initial_size-len(G)
		controlled_nodes_slice, matching = get_controlled_nodes(G)
		(source, external, internal), types = control_profile_(G, controlled_nodes_slice, matching)
		if not remove_disconnected_nodes:
			controlled_nodes_binary = np.ones((len(G), 1))
			for nidx in controlled_nodes_slice:
				if types[nidx] == TYPE_SOURCE:
					controlled_nodes_binary[nidx] = 0.9
				elif types[nidx] == TYPE_EXTERNAL_DILATION:
					controlled_nodes_binary[nidx] = 0.3
				else:
					controlled_nodes_binary[nidx] = 0
		for j in range(0, len(controlled_nodes_binary)):
			controlled_nodes[j][i] = controlled_nodes_binary[j]
		m = len(controlled_nodes_slice)
		controls[i] = (m, source, external, internal)
		if i > 0:
			diff_controls[i-1] = m-controls[i-1][0]
	print 'Ran in %f seconds' % (time.time()-total_running_time)
	return controls, diff_controls, weight_cuts, node_removed, controlled_nodes

def print_control_graphs(controls, rcontrols, filename):
	N = len(controls[:,0])
	ind = np.arange(N)
	width = 1

	fig = plt.figure()
	plt.subplot(2, 2, 1)
	p1 = plt.bar(ind, controls[:,1], width, color='y', edgecolor="none", linewidth=0)
	p2 = plt.bar(ind, controls[:,2], width, color='b', bottom=controls[:,1], edgecolor="none", linewidth=0)
	p3 = plt.bar(ind, controls[:,3], width, color='k', bottom=controls[:,1]+controls[:,2], edgecolor="none", linewidth=0)
	plt.title('Weight cuts')
	plt.xlabel('slices')
	plt.ylabel('fraction of controls')
	plt.legend( (p1[0], p2[0], p3[0]), ('Source controls', 'External dilation', 'Internal dilation'), bbox_to_anchor=(0., 1.08, 1., .102), loc=3, ncol=1, borderaxespad=0. )
	# Control profiles on random cuts
	plt.subplot(2, 2, 2)
	plt.bar(ind, rcontrols[:,1], width, color='y', edgecolor="none", linewidth=0)
	plt.bar(ind, rcontrols[:,2], width, color='b', bottom=rcontrols[:,1], edgecolor="none", linewidth=0)
	plt.bar(ind, rcontrols[:,3], width, color='k', bottom=rcontrols[:,1]+rcontrols[:,2], edgecolor="none", linewidth=0)
	plt.title('Random cuts')
	plt.xlabel('slices')
	plt.ylabel('fraction of controls')
	# Absolute number of controls 
	plt.subplot(2, 2, 3)
	plt.bar(ind, controls[:,0]*controls[:,1], width, color='y', edgecolor="none", linewidth=0)
	plt.bar(ind, controls[:,0]*controls[:,2], width, color='b', bottom=controls[:,0]*controls[:,1], edgecolor="none", linewidth=0)
	plt.bar(ind, controls[:,0]*controls[:,3], width, color='k', bottom=controls[:,0]*controls[:,1]+controls[:,0]*controls[:,2], edgecolor="none", linewidth=0)
	plt.title('Weight cuts')
	plt.xlabel('slices')
	plt.ylabel('Number of controls')
	plt.subplot(2, 2, 4)
	plt.bar(ind, rcontrols[:,0]*rcontrols[:,1], width, color='y', edgecolor="none", linewidth=0)
	plt.bar(ind, rcontrols[:,0]*rcontrols[:,2], width, color='b', bottom=rcontrols[:,0]*rcontrols[:,1], edgecolor="none", linewidth=0)
	plt.bar(ind, rcontrols[:,0]*rcontrols[:,3], width, color='k', bottom=rcontrols[:,0]*rcontrols[:,1]+rcontrols[:,0]*rcontrols[:,2], edgecolor="none", linewidth=0)
	plt.title('Random cuts')
	plt.xlabel('slices')
	plt.ylabel('Number of controls')
	fig.savefig('graphs/'+filename+'-controls.png')
	plt.show()

def print_graphs(controls, rcontrols, diff_controls, weight_cuts, node_removed, rnode_removed, filename):
	N = len(controls[:,0])
	ind = np.arange(N)
	width = 1
	# Number of controls, weight cuts vs random cuts
	fig = plt.figure()
	plt.subplot(3, 1, 1)
	plt.title("Minimum number of controls")
	plt.plot(range(0,len(controls[:,0])), controls[:,0], 'b-', label="Weight cuts")
	plt.plot(range(0,len(controls[:,0])), rcontrols[:,0], 'g-', label="Random cuts")
	plt.ylabel('controls')
	plt.legend(loc=2)
	# Weight of the heaviest edge in the cut
	plt.subplot(3, 1, 2)
	plt.title("Weight of the heaviest edge in the cut")
	plt.plot(range(1,len(weight_cuts)+1), weight_cuts, 'b-', label="Weight of edge cut")
	plt.ylabel('weight')
	plt.legend(loc=2)
	# Differential of controls between each slice
	plt.subplot(3, 1, 3)
	plt.title("Differential of controls between each slice")
	plt.plot(range(1,len(diff_controls)+1), diff_controls, 'b-')
	plt.axis([0, len(diff_controls)+1, -0.1, max(diff_controls)+0.5])
	plt.xlabel('slices')
	fig.savefig('graphs/'+filename+'-cuts.png')
	plt.show()

	# Control profiles on weight cuts
	fig = plt.figure()
	plt.subplot(2, 2, 1)
	p1 = plt.bar(ind, controls[:,1], width, color='y', edgecolor="none", linewidth=0)
	p2 = plt.bar(ind, controls[:,2], width, color='b', bottom=controls[:,1], edgecolor="none", linewidth=0)
	p3 = plt.bar(ind, controls[:,3], width, color='k', bottom=controls[:,1]+controls[:,2], edgecolor="none", linewidth=0)
	plt.title('Weight cuts')
	plt.xlabel('slices')
	plt.ylabel('fraction of controls')
	plt.legend( (p1[0], p2[0], p3[0]), ('Source controls', 'External dilation', 'Internal dilation'), bbox_to_anchor=(0., 1.08, 1., .102), loc=3, ncol=1, borderaxespad=0. )
	# Control profiles on random cuts
	plt.subplot(2, 2, 2)
	plt.bar(ind, rcontrols[:,1], width, color='y', edgecolor="none", linewidth=0)
	plt.bar(ind, rcontrols[:,2], width, color='b', bottom=rcontrols[:,1], edgecolor="none", linewidth=0)
	plt.bar(ind, rcontrols[:,3], width, color='k', bottom=rcontrols[:,1]+rcontrols[:,2], edgecolor="none", linewidth=0)
	plt.title('Random cuts')
	plt.xlabel('slices')
	plt.ylabel('fraction of controls')
	# Absolute number of controls 
	plt.subplot(2, 2, 3)
	plt.bar(ind, controls[:,0]*controls[:,1], width, color='y', edgecolor="none", linewidth=0)
	plt.bar(ind, controls[:,0]*controls[:,2], width, color='b', bottom=controls[:,0]*controls[:,1], edgecolor="none", linewidth=0)
	plt.bar(ind, controls[:,0]*controls[:,3], width, color='k', bottom=controls[:,0]*controls[:,1]+controls[:,0]*controls[:,2], edgecolor="none", linewidth=0)
	plt.title('Weight cuts')
	plt.xlabel('slices')
	plt.ylabel('Number of controls')
	plt.subplot(2, 2, 4)
	plt.bar(ind, rcontrols[:,0]*rcontrols[:,1], width, color='y', edgecolor="none", linewidth=0)
	plt.bar(ind, rcontrols[:,0]*rcontrols[:,2], width, color='b', bottom=rcontrols[:,0]*rcontrols[:,1], edgecolor="none", linewidth=0)
	plt.bar(ind, rcontrols[:,0]*rcontrols[:,3], width, color='k', bottom=rcontrols[:,0]*rcontrols[:,1]+rcontrols[:,0]*rcontrols[:,2], edgecolor="none", linewidth=0)
	plt.title('Random cuts')
	plt.xlabel('slices')
	plt.ylabel('Number of controls')
	fig.savefig('graphs/'+filename+'-controls.png')
	plt.show()

def avg_degrees_(G, nidx):
	total = 0.0
	total_in = 0.0
	total_out = 0.0
	avg = 0.0
	avg_in = 0.0
	avg_out = 0.0
	for eidx, w in G.edges_iter_(nidx, weight=True):
		total += w
	for eidx, w in G.in_edges_iter_(nidx, weight=True):
		total_in += w
	for eidx, w in G.out_edges_iter_(nidx, weight=True):
		total_out += w
	if G.degree_(nidx) > 0:
		avg = total/G.degree_(nidx)
	if G.in_degree_(nidx) > 0:
		avg_in = total_in/G.in_degree_(nidx)
	if G.out_degree_(nidx) > 0:
		avg_out = total_out/G.out_degree_(nidx)
	return avg, avg_in, avg_out

def weight_means(G):
	total = np.zeros(len(G))
	total_in = np.zeros(len(G))
	total_out = np.zeros(len(G))
	for nidx in G.nodes_iter_():
		(avg, avg_in, avg_out) = avg_degrees_(G, nidx)
		total[nidx] = avg
		total_in[nidx] = avg_in
		total_out[nidx] = avg_out
	return total, total_in, total_out

def scatter_plot(x, y):
	fig = plt.figure()
	plt.scatter(x, y)
	plt.show()

def show_weights(G):
	(total, total_in, total_out) = weight_means(G)
	(deg, in_deg, out_deg) = degree_distribution(G)
	fig = plt.figure()
	plt.subplot(3,1,1)
	plt.scatter(deg, total)
	plt.subplot(3,1,2)
	plt.scatter(in_deg, total_in)
	plt.subplot(3,1,3)
	plt.scatter(out_deg, total_out)
	plt.show()	

def correlations(G):
	(total, total_in, total_out) = weight_means(G)
	(deg, in_deg, out_deg) = degree_distribution(G)
	return np.corrcoef(total, deg)[0,1], np.corrcoef(total_in, in_deg)[0,1], np.corrcoef(total_out, out_deg)[0,1]
		
def tag_nodes(G):
	controlled_nodes_slice, matching = get_controlled_nodes(G)
	(source, external, internal), types = control_profile_(G, controlled_nodes_slice, matching)
	return types

def filename(string):
	return string.split("/")[-1].split('.')[0]


	
