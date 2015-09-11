import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
import numpy as np
import zen
import random
import os
import time
import pandas as pd
import ggplot as gg
import powerlaw
import ternary

def load_network(networkname):
	G = zen.io.edgelist.read(networkname, directed=True, weighted=True)
	return G

def random_barabasi():
    G = zen.generating.barabasi_albert(100, 10, directed=True)
    for eidx in G.edges_iter_():
        (src, tgt) = G.endpoints_(eidx)
        # Increasing weights
        G.set_weight_(eidx, random.randint(1,100)*(G.degree_(src)+G.degree_(tgt)))
        # Decreasing weights
        # G.set_weight_(eidx, random.randint(1,100)*(G.degree_(src)+G.degree_(tgt)))
        # Random weights
        # G.set_weight_(eidx, random.randint(1,100))
    return G
    
def random_erdos():
    G = zen.generating.erdos_renyi(100, 0.02, directed=True)
    for eidx in G.edges_iter_():
        (src, tgt) = G.endpoints_(eidx)
        # Increasing weights
        # G.set_weight_(eidx, random.randint(1,100)*(G.degree_(src)+G.degree_(tgt)))
        # Decreasing weights
        G.set_weight_(eidx, random.randint(1,100)/(G.degree_(src)+G.degree_(tgt)))
        # Random weights
        # G.set_weight_(eidx, random.randint(1,100))
    return G

def filename(string):
	return string.split("/")[-1].split('.')[0]

def degree_distribution(G):
	in_degrees = np.array([G.in_degree_(nidx) for nidx in G.nodes_iter_()])
	out_degrees = np.array([G.out_degree_(nidx) for nidx in G.nodes_iter_()])
	degrees = np.array([G.degree_(nidx) for nidx in G.nodes_iter_()])
	return (in_degrees, out_degrees, degrees)

def plot_deg_distrib(G):
	(in_deg, out_deg, deg) = degree_distribution(G)
	in_deg_series = pd.Series(in_deg)
	out_deg_series = pd.Series(out_deg)
	in_out = { 'in_deg': in_deg_series, 'out_deg': out_deg_series }
	df = pd.DataFrame(in_out)
	df = pd.melt(df)
	p = gg.ggplot(gg.aes(x='value', color='variable', fill='variable'), data=df2) + gg.geom_histogram(alpha=0.6, binwidth=1)
	print p

def avg_weights(G):
    (indeg, outdeg, deg) = degree_distribution(G)
    maxin = max(indeg)+1
    maxout = max(outdeg)+1
    total_max = max(deg)+1
    number_in = np.zeros((maxin, maxout)) # (k, l) is the number of in-neighbors for nodes of type (k, l)
    number_out = np.zeros((maxin, maxout)) # (k, l) is the number of out-neighbors for nodes of type (k, l)
    total_number = np.zeros(total_max) # (k) is the number of neighbors for nodes of degree k
    weight_in = np.zeros((maxin, maxout)) # (k, l) is the total weight of in-edges for nodes of type (k, l)
    weight_out = np.zeros((maxin, maxout)) # (k, l) is the total weight of out-edges for nodes of type (k, l)
    total_weight = np.zeros(total_max) # (k) is the total weight of edges for nodes of degree k
    avg_in = np.zeros((maxin, maxout)) # (k, l) is the average weight of in-edges for nodes of type (k, l)
    avg_out = np.zeros((maxin, maxout)) # (k, l) is the average weight of out-edges for nodes of type (k, l)
    avg = np.zeros((maxin, maxout)) # (k, l) is the average weight of edges for nodes of type (k, l)
    total_avg_in = np.zeros(maxin) # (k) is the average weight of edges for nodes of in-degree k
    total_avg_out = np.zeros(maxout) # (k) is the average weight of edges for nodes of out-degree k
    total_avg = np.zeros(total_max) # (k) is the average weight of edges for nodes of degree k

    for nidx in G.nodes_iter_():
    	node_out_deg = G.out_degree_(nidx)
    	node_in_deg = G.in_degree_(nidx)
    	node_deg = G.degree_(nidx)
    	number_in[node_in_deg, node_out_deg] += node_in_deg
    	number_out[node_in_deg, node_out_deg] += node_out_deg
    	total_number[node_deg] += node_deg
    	for eidx, w in G.in_edges_iter_(nidx, weight=True): # edges for which nidx is a target
    		weight_in[node_in_deg, node_out_deg] += w
    		total_weight[node_deg] += w
            
    	for eidx, w in G.out_edges_iter_(nidx, weight=True):
    		weight_out[node_in_deg, node_out_deg] += w
    		total_weight[node_deg] += w
            
    for k in range(0, maxin):
    	for l in range(0, maxout):
    		if number_in[k, l] > 0:
    			avg_in[k,l] = weight_in[k,l]/number_in[k,l]
    		if number_out[k, l] > 0:
    			avg_out[k,l] = weight_out[k,l]/number_out[k,l]
    		if number_in[k,l] > 0 or number_out[k,l] > 0:
    			avg[k,l] = (weight_in[k,l]+weight_out[k,l])/(number_in[k,l]+number_out[k,l])

    for k in range(0, total_max):
    	if total_number[k] > 0:
    		total_avg[k] = total_weight[k] / total_number[k]       

    return (avg_in, avg_out, avg, total_avg_in, total_avg_out, total_avg)

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
    # Runs reducer <runs> times and averages the results
    # to erase randomness in the choice of links
	iter_number = H.size()/step + 1
	avg_controls = np.zeros((iter_number, 4))
	for i in range(0, runs):
		print 'Run %d' % i
		[controls, diff_controls, weight_cuts, node_removed, controlled_nodes] = reducer(H, step, rand, asc)
		avg_controls = avg_controls+controls/runs
	if filename != "":
		np.savetxt("controls/"+filename+".txt", avg_controls)
	return avg_controls

def reducer(H, step, rand=False, asc=True):
    # Removes edges <step> at a time and keeps track of the evolution of control profile
    # rand=False -> weight cuts
    # rand=True -> random cuts
    # asc=True -> Lighter edges are cut first
	total_running_time = time.time()
	G = H.copy()
	w = [w for u,v,w in G.edges_iter(weight=True)]
	iter_number = G.size()/step + 1
	controls = np.zeros((iter_number, 4))	
	diff_controls = np.zeros(iter_number-1)
	weight_cuts = np.zeros(iter_number-1)
	controlled_nodes = np.ones((len(G), iter_number))
	remove_disconnected_nodes = False # remove disconnected nodes
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
    # controls holds the control profiles at each slice
    # diff_controls holds the added number of controls between each slice
    # weight_cuts holds the max weight cut during each slice
    # node_removed has number of nodes removed until slice
    # controlled_nodes gives the indices of controlled nodes at each slice
	return controls, diff_controls, weight_cuts, node_removed, controlled_nodes

def correlations(G):
    (avg_in, avg_out, avg, total_avg_in, total_avg_out, total_avg) = avg_weights(G)
    (in_deg, out_deg, deg) = degree_distribution(G)
    nz = total_avg[np.nonzero(total_avg)]
    r = np.array(range(0, len(total_avg)))
    r = r[np.nonzero(total_avg)]
    return np.corrcoef(nz, r)[0,1]
		
def plot_link_weights(G, filename):
	print "Filename = %s" % filename
	w = np.array([w for u,v,w in G.edges_iter(weight=True)])
	fit = powerlaw.Fit(w)
	R, p = fit.distribution_compare('power_law', 'exponential', normalized_ratio=True)
	title = "Exponent: %f" % fit.power_law.alpha
	print "(R, p) = (%f, %f)" % (R, p)
	df = pd.DataFrame({ 'weight': pd.Series(w) })
	p = gg.ggplot(gg.aes(x='weight'), data=df) + gg.geom_histogram(fill='blue', alpha=0.6) + gg.scale_y_log10() + gg.ggtitle(title)
	gg.ggsave(filename="graphs/weight_distrib/"+filename+"-weight_distrib.png", plot=p)
	
def get_exponent(G):
    w = np.array([w for u,v,w in G.edges_iter(weight=True)])
    fit = powerlaw.Fit(w)
    R, p = fit.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    return fit.power_law.alpha

def draw_ternary_plot(controls, rcontrols, filename):
	figure, tax = ternary.figure(scale=1.0)
	tax.boundary(color='black')
	tax.plot(controls, linewidth=2.0, label="Weight cuts", color="blue")
	tax.plot(rcontrols, linewidth=2.0, label="Random cuts", color="green")
	tax.left_axis_label("No source control")
	tax.right_axis_label("No internal dilation control")
	tax.bottom_axis_label("No external dilation control")
	tax.legend()
	#tax.show()
	figure.savefig('graphs/ternary/'+filename+'-ternary.png')

def plot_min_number_controls(controls, rcontrols, filename):
	fig = plt.figure()
	plt.plot(range(0, len(controls)), controls, 'b-', label='Weight cuts', linewidth=2.0)
	plt.plot(range(0, len(rcontrols)), rcontrols, 'g-', label='Random cuts', linewidth=2.0)
	plt.legend(loc=2)
	plt.title(filename)
	plt.ylabel('Number of controls')
	plt.xlabel('Slices')
#	plt.show()
	fig.savefig('graphs/number_controls/'+filename+'-controls.png')




