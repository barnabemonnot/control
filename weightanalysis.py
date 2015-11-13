import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import zen
import random
import os
import time
import pandas as pd
import ggplot as gg
from powerlaw import *
import ternary
from math import *
from scipy.stats import *
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go

def swap_cols(arr, frm, to):
    arr[:,[frm, to]] = arr[:,[to, frm]]

def filename(string):
	return string.split("/")[-1].split('.')[0]

def load_network(networkname):
	G = zen.io.edgelist.read(networkname, directed=True, weighted=True)
	return G

def get_weights(G):
    w = np.array([w for u,v,w in G.edges_iter(weight=True)])
    return np.sort(w)

def random_barabasi():
    G = zen.generating.barabasi_albert(1000, 2, directed=True)
    return G

def random_erdos():
    b = True
    n = 1000
    k = 2.0
    G = []
    while b:
        G = zen.generating.erdos_renyi(n, k/(n-1), directed=True)
        print "Generating one"
        print 2*len(G)
        print G.size()
        if abs(G.size()-2*n) <= 0.0001*k*n:
            b = False
    return G

def assign_weights(G, mode):
    for eidx in G.edges_iter_():
        (src, tgt) = G.endpoints_(eidx)
        if mode == "inc":
            # Increasing weights
            G.set_weight_(eidx, float(random.randint(1,100))*(G.degree_(src)+G.degree_(tgt)))
        elif mode == "inv":
            # Decreasing weights
            G.set_weight_(eidx, float(random.randint(1,100))/(G.degree_(src)+G.degree_(tgt)))
        else:
            # Random weights
            G.set_weight_(eidx, float(random.randint(1,100)))

def assign_weights_in(G, mode):
    for eidx in G.edges_iter_():
        (src, tgt) = G.endpoints_(eidx)
        if mode == "inc":
            G.set_weight_(eidx, float(random.randint(1,100))*G.in_degree_(tgt))
        elif mode == "inv":
            G.set_weight_(eidx, float(random.randint(1,100))/G.in_degree_(tgt))

def assign_weights_out(G, mode):
    for eidx in G.edges_iter_():
        (src, tgt) = G.endpoints_(eidx)
        if mode == "inc":
            G.set_weight_(eidx, float(random.randint(1,100))*G.out_degree_(src))
        elif mode == "inv":
            G.set_weight_(eidx, float(random.randint(1,100))/G.out_degree_(src))

def assign_weights_in_out(G, mode):
    for eidx in G.edges_iter_():
        (src, tgt) = G.endpoints_(eidx)
        if mode == "inv":
            G.set_weight_(eidx, float(random.randint(1,100))/G.out_degree_(src)/G.in_degree_(tgt))
        elif mode == "inc":
            G.set_weight_(eidx, float(random.randint(1,100))*G.out_degree_(src)*G.in_degree_(tgt))

def assign_weights_in_out_rev(G, mode):
    for eidx in G.edges_iter_():
        (src, tgt) = G.endpoints_(eidx)
        if mode == "inv":
            G.set_weight_(eidx, float(random.randint(1,100))*G.out_degree_(src)/G.in_degree_(tgt))
        elif mode == "inc":
            G.set_weight_(eidx, float(random.randint(1,100))/G.out_degree_(src)*G.in_degree_(tgt))

def degree_distribution(G):
	in_degrees = np.array([G.in_degree_(nidx) for nidx in G.nodes_iter_()])
	out_degrees = np.array([G.out_degree_(nidx) for nidx in G.nodes_iter_()])
	degrees = np.array([G.degree_(nidx) for nidx in G.nodes_iter_()])
	return (in_degrees, out_degrees, degrees)

def avg_control_distance(controls, rcontrols):
    return np.sum((controls-rcontrols)/np.max(controls))/len(controls)

def avg_weights_degree(G):
    avg_weight = []
    avg_weight_in = []
    avg_weight_out = []
    deg = []
    in_deg = []
    out_deg = []
    for i in range(0, len(G)):
        if G.degree_(i) > 0:
            avg_weight += [sum([w for eidx,w in G.edges_iter_(i, weight=True)])/float(G.degree_(i))]
            deg += [G.degree_(i)]
        if G.in_degree_(i) > 0:
            avg_weight_in += [sum([w for eidx,w in G.in_edges_iter_(i, weight=True)])/float(G.in_degree_(i))]
            in_deg += [G.in_degree_(i)]
        if G.out_degree_(i) > 0:
            avg_weight_out += [sum([w for eidx,w in G.out_edges_iter_(i, weight=True)])/float(G.out_degree_(i))]
            out_deg += [G.out_degree_(i)]
    return (np.array(avg_weight_in), np.array(in_deg), np.array(avg_weight_out), np.array(out_deg), np.array(avg_weight), np.array(deg))

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

def avg_reducer(H, step, runs=100, rand=False, asc=True, filename="", log=False):
    # Runs reducer <runs> times and averages the results
    # to erase randomness in the choice of links
    iter_number = H.size()/step + 1
    avg_controls = np.zeros((iter_number, 4))
    total_weight_cuts = np.zeros(iter_number-1)
    for i in range(0, runs):
        if log:
            print 'Run %d' % i
        [controls, diff_controls, weight_cuts, node_removed, controlled_nodes] = reducer(H, step, rand=rand, asc=asc, log=log)
        avg_controls = avg_controls+controls/runs
        total_weight_cuts = total_weight_cuts+weight_cuts/runs
        if filename != "":
            np.savetxt("controls/"+filename+".txt", avg_controls)
    return (avg_controls, total_weight_cuts)

def reducer(H, step, rand=False, asc=True, log=False):
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
        if i % 100 == 0 and log:
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
    if log:
        print 'Ran in %f seconds' % (time.time()-total_running_time)
    # controls holds the control profiles at each slice
    # diff_controls holds the added number of controls between each slice
    # weight_cuts holds the max weight cut during each slice
    # node_removed has number of nodes removed until slice
    # controlled_nodes gives the indices of controlled nodes at each slice
    return controls, diff_controls, weight_cuts, node_removed, controlled_nodes

def plot_deg_distrib(G):
	(in_deg, out_deg, deg) = degree_distribution(G)
	in_deg_series = pd.Series(in_deg)
	out_deg_series = pd.Series(out_deg)
	in_out = { 'in_deg': in_deg_series, 'out_deg': out_deg_series }
	df = pd.DataFrame(in_out)
	df = pd.melt(df)
	p = gg.ggplot(gg.aes(x='value', color='variable', fill='variable'), data=df2) + gg.geom_histogram(alpha=0.6, binwidth=1)
	print p

def spearman_correlation(G):
    (avg_in, deg_in, avg_out, deg_out, avg, deg) = avg_weights_degree(G)
    return (spearmanr(np.array([avg_in,deg_in]), axis=1)[0], spearmanr(np.array([avg_out,deg_out]), axis=1)[0], spearmanr(np.array([avg,deg]), axis=1)[0])

def empirical_correlation(G):
    (avg_in, deg_in, avg_out, deg_out, avg, deg) = avg_weights_degree(G)
    return (np.corrcoef(avg_in, deg_in)[0,1], np.corrcoef(avg_out, deg_out)[0,1], np.corrcoef(avg, deg)[0,1])

def get_exponent(G):
    w = np.array([w for u,v,w in G.edges_iter(weight=True)])
    fit = Fit(w)
    R, p = fit.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    #print "R = %f; p = %f" % (R, p)
    return fit.power_law.alpha

def plot_hist_ly(G, filename):
    w = np.array([w for u,v,w in G.edges_iter(weight=True)])
    trace = go.Histogram(
        x=w,
        histnorm=1
    )
    data = [trace]
    layout = go.Layout(
        title='<b>Weight distribution for ' + filename + '</b>',
        font=dict(size=18),
        xaxis=dict(
            title='Weight',
            type='log'
        ),
        yaxis=dict(
            title='Frequency',
            type='log',
            autorange=True
        ),
        bargap=0.25,
        bargroupgap=0.3
    )
    fig = go.Figure(data=data, layout=layout)
    plot_url = py.plot(fig, filename='basic-line')
    print plot_url
	
def plot_link_weights(G, filename):
    print "Filename = %s" % filename
    w = np.array([w for u,v,w in G.edges_iter(weight=True)])
    print w
    fit = Fit(w)
    R, p = fit.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    title = "Exponent: %f" % fit.power_law.alpha
    print "(R, p) = (%f, %f)" % (R, p)
    df = pd.DataFrame({ 'weight': pd.Series(np.log10(w)) })
    p = gg.ggplot(gg.aes(x='weight'), data=df) + gg.geom_histogram(fill='blue', alpha=0.6) + gg.scale_y_log10() + gg.ggtitle(title)
    gg.ggsave(filename="graphs/weight_distrib/"+filename+"-weight_distrib.png", plot=p)

def scatter(x, y, filename=""):
    df = pd.DataFrame({ 'x': pd.Series(x), 'y': pd.Series(y) })
    p = gg.ggplot(gg.aes(x='x', y='y'), data=df) + gg.geom_point()
    if filename == "":
        print p
    else:
        gg.ggsave(filename="graphs/scatter/"+filename+".png", plot=p)
    
def scatter_weights_degree(G, filename=""):
    (avg, deg) = avg_weights_degree(G)
    scatter(avg, deg, filename)
    
def plot_min_ly(controls, rcontrols, weights, filename):
    trace1 = go.Scatter(
        x = np.array(range(0, len(controls)))/float(len(controls)-1)*100.0,
        y = controls/np.max(controls),
        name = "Weight cuts",
        line = dict(color='#45A9DE',width=4)
    )
    
    trace2 = go.Scatter(
        x = np.array(range(0, len(rcontrols)))/float(len(rcontrols)-1)*100.0,
        y = rcontrols/np.max(rcontrols),
        name = "Random cuts",
        line = dict(dash='dashdot',color='#BFCE80',width=4)
    )
    
    trace3 = go.Scatter(
        x = np.array(range(0, len(rcontrols)))/float(len(rcontrols)-1)*100.0,
        y = weights,
        line = dict(dash='dot', color='#45A9DE',width=4),
        showlegend=False
    )
    
    fig = tools.make_subplots(rows=3, cols=1, specs=[[{'rowspan': 2}],[None],[{}]])
    fig['layout'].update(title='Minimum number of controls for ' + filename, font=dict(size=18), height=600, width=800)
    fig['layout']['yaxis1'].update(title='Fraction of driver nodes', range=[0,1], titlefont=dict(size=16))
    fig['layout']['yaxis2'].update(title='Sliced edge weight', titlefont=dict(size=16))
    fig['layout']['xaxis2'].update(title='Percentage of slices', titlefont=dict(size=16))
    
    fig.append_trace(trace1,1,1)
    fig.append_trace(trace2,1,1)
    fig.append_trace(trace3,3,1)

    plot_url = py.plot(fig, filename='basic-line')
    print plot_url

def plot_min_number_controls(controls, rcontrols, filename):
    fig = plt.figure()
    plt.plot(np.array(range(0, len(controls)))/float(len(controls)-1)*100.0, controls/np.max(controls), 'b-', label='Weight cuts', linewidth=2.0)
    plt.plot(np.array(range(0, len(rcontrols)))/float(len(controls)-1)*100.0, rcontrols/np.max(rcontrols), 'g--', label='Random cuts', linewidth=2.0)
    plt.legend(loc=2)
    plt.title(filename)
    plt.axis([0, 100, 0, 1])
    plt.ylabel('Fraction of controlled nodes')
    plt.xlabel('Percentage of edge cuts')
    plt.grid()
    plt.show()
    fig.savefig('graphs/number_controls/'+filename+'-controls.png')

def draw_ternary_plot(controls, rcontrols, filename):
    m = np.shape(controls)[0]
    n = 20
    ten_per = int(floor(m/5.0))
    swap_cols(controls, 2, 0)
    swap_cols(rcontrols, 2, 0)
    swap_cols(controls, 0, 1) 
    swap_cols(rcontrols, 0, 1)
    first_ten_per = controls[range(0, ten_per),:]
    rfirst_ten_per = rcontrols[range(0, ten_per),:]
    rest = controls[range(ten_per, m),:]
    rrest = rcontrols[range(ten_per, m),:]
    figure, tax = ternary.figure(scale=1.0)
    tax.boundary(color='black')
    tax.plot(first_ten_per, linewidth=3, label="Weight cuts", color="#00d3e4")
    tax.plot(rfirst_ten_per, linewidth=3, label="Random cuts", color="#75b765")
    tax.plot(rest, linewidth=1, label="Weight cuts", color="#007e88")
    tax.plot(rrest, linewidth=1, label="Random cuts", color="#466d3c")
    tax.left_axis_label("No external dilation control")
    tax.right_axis_label("No source control")
    tax.bottom_axis_label("No internal dilation control")
    tax.legend()
    #tax.show()
    figure.savefig('graphs/ternary/celegans/'+filename+'-ternary.png')



