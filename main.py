import weightanalysis as wa
import argparse
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser()
parser.add_argument('--n', '--network', type=str, default='', help='network', metavar='n')
args = parser.parse_args()

G = wa.load_network(args.n)
step = 20
#(deg, in_deg, out_deg) = wa.degree_distribution(G)
#print "Will plot degree distribution"
#wa.plot_deg_distrib(G)
#print "Compute weight cuts"
#(con, diff, weight, node, controlled) = wa.reducer(G, step, asc=False)
#print "Compute random cuts"
#(rcon, rdiff, rweight, rnode, rcontrolled) = wa.reducer(G, step, rand=True, asc=False)
#print "Print graphs"
#filename = wa.filename(args.n)
#wa.print_graphs(con, rcon, diff, weight, node, rnode, filename)
#(total, total_in, total_out) = wa.weight_means(G)
#wa.scatter_plot(deg, total)
#wa.show_weights(G)
#(corr_deg, corr_in, corr_out) = wa.correlations(G)
#print "Network %s, correlations deg = %f, in = %f, out = %f" % (args.n, corr_deg, corr_in, corr_out)
print 'Avg_reducer weight'
filename = wa.filename(args.n)
avg = wa.avg_reducer(G, step, runs=20, filename=filename)
print 'Avg_reducer random'
r_avg = wa.avg_reducer(G, step, rand=True, runs=20, filename=filename+"rand")
wa.print_control_graphs(avg, r_avg, wa.filename(args.n)+'avg')
