import weightanalysis as wa
import argparse
import numpy as np
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--n', '--network', type=str, default='', help='network', metavar='n')
args = parser.parse_args()

# filename = wa.filename(args.n)
# print filename
# G = wa.load_network(args.n)
filename="erinv"
G = wa.random_erdos()
avg = wa.avg_reducer(G, 1, filename=filename)
r_avg = wa.avg_reducer(G, 1, rand=True, filename=filename+"-rand")
# desc = {
#     "intl": "International airport connections",
#     "us": "US airport connections",
#     "celegans": "Neural network of \\textit{C.Elegans} bacteria",
#     "Chesapeake": "Foodwebs",
#     "ChesLower": "Foodwebs",
#     "ChesMiddle": "Foodwebs",
#     "ChesUpper": "Foodwebs",
#     "CrystalC": "Foodwebs",
#     "CrystalD": "Foodwebs",
#     "Everglades": "Foodwebs",
#     "Florida": "Foodwebs",
#     "Maspalomas": "Foodwebs",
#     "Michigan": "Foodwebs",
#     "Mondego": "Foodwebs",
#     "Narragan": "Foodwebs",
#     "Rhode": "Foodwebs",
#     "StMarks": "Foodwebs",
#     "days1-12": "Mail logs (all days)",
#     "days1-4": "Mail logs",
#     "days5-8": "Mail logs",
#     "days9-12": "Mail logs",
#     "mention_network": "Twitter network of mentions",
#     "one_mode_char": "ucirvine",
#     "one_mode_message": "ucirvine"
# }
#
# # Name / Description / #Node / #Edges (/ gamma) / rho /
# print "\hline"
# print "%s & %s & %d & %d & %f & %f \\\\" % (filename, desc[filename], len(G), G.size(), wa.get_exponent(G), wa.correlations(G))
# step = 1

#print 'Avg_reducer weight'
#avg = wa.avg_reducer(G, step, runs=20, filename=filename)
#print 'Avg_reducer random'
#r_avg = wa.avg_reducer(G, step, rand=True, runs=20, filename=filename+"-rand")

# controls = np.loadtxt("controls/"+filename+".txt")
# profile = controls[:,range(1,4)]
# rcontrols = np.loadtxt("controls/"+filename+"-rand.txt")
# rprofile = rcontrols[:,range(1,4)]
#wa.draw_ternary_plot(profile, rprofile, filename)
wa.plot_min_number_controls(avg[:,0], r_avg[:,0], filename)
#wa.gg_print_control(controls)
#(avg_in, avg_out, avg, total_in, total_out, total) = wa.avg_weights(G)
#for nidx in G.nodes_iter_():
#	if G.out_degree_(nidx) == 0:
#		#print "I am a sink with %d neighbors" % G.in_degree_(nidx)
#		for eidx, w in G.edges_iter_(nidx, weight=True):
#			#print "I have weights %f" % w
#print "Compute weight cuts"
#(con, diff, weight, node, controlled) = wa.reducer(G, step, asc=False)
#print "Compute random cuts"
#(rcon, rdiff, rweight, rnode, rcontrolled) = wa.reducer(G, step, rand=True, asc=False)

#(deg, in_deg, out_deg) = wa.degree_distribution(G)
#print "Will plot degree distribution"
#wa.plot_deg_distrib(G)
#wa.plot_link_weights(G, filename)
#print "Print graphs"
#filename = wa.filename(args.n)
#wa.print_graphs(con, rcon, diff, weight, node, rnode, filename)
#(total, total_in, total_out) = wa.weight_means(G)
#wa.scatter_plot(deg, total)
#wa.show_weights(G)
#(corr_deg, corr_in, corr_out) = wa.correlations(G)
#print "Network %s, correlations deg = %f, in = %f, out = %f" % (args.n, corr_deg, corr_in, corr_out)
#wa.print_control_graphs(avg, r_avg, wa.filename(args.n)+'avg')
