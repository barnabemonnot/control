import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
import numpy as np
import ternary
import weightanalysis as wa
import pandas as pd
import ggplot as gg
from powerlaw import *
from scipy.stats import *
from math import *
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go

def swap_cols(arr, frm, to):
    arr[:,[frm, to]] = arr[:,[to, frm]]

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

def plot_deg_distrib(G):
	(in_deg, out_deg, deg) = wa.degree_distribution(G)
	in_deg_series = pd.Series(in_deg)
	out_deg_series = pd.Series(out_deg)
	in_out = { 'in_deg': in_deg_series, 'out_deg': out_deg_series }
	df = pd.DataFrame(in_out)
	df = pd.melt(df)
	p = gg.ggplot(gg.aes(x='value', color='variable', fill='variable'), data=df2) + gg.geom_histogram(alpha=0.6, binwidth=1)
	print p

def plot_hist_ly(G, filename):
    w = np.array([w for u,v,w in G.edges_iter(weight=True)])
    trace = go.Histogram(
        x=np.log10(w),
        histnorm=1,
        marker=dict(
            color='#45A9DE'
        )
    )
    data = [trace]
    layout = go.Layout(
        font=dict(size=18),
        xaxis=dict(
            title='Weight',
            titlefont=dict(size=28),
            tickfont=dict(size=24)
        ),
        yaxis=dict(
            title='Frequency',
            type='log',
            autorange=True,
            titlefont=dict(size=28),
            tickfont=dict(size=24)
        ),
        bargap=0.01,
        bargroupgap=0.01,
    )
    fig = go.Figure(data=data, layout=layout)
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
    figure.savefig('graphs/ternary/'+filename+'-ternary.png')

def plot_min_ly(controls, rcontrols, weights, filename):
    trace1 = go.Scatter(
        x = np.array(range(0, len(controls)))/float(len(controls)-1)*100.0,
        y = controls/np.max(controls),
        name = "Weight cuts",
        line = dict(color='#45A9DE',width=4),
        showlegend=False
    )

    trace2 = go.Scatter(
        x = np.array(range(0, len(rcontrols)))/float(len(rcontrols)-1)*100.0,
        y = rcontrols/np.max(rcontrols),
        name = "Random cuts",
        line = dict(dash='dashdot',color='#BFCE80',width=4),
        showlegend=False
    )

    trace3 = go.Scatter(
        x = np.array(range(0, len(weights)))/float(len(weights)-1)*100.0,
        y = weights,
        line = dict(dash='dot', color='#45A9DE',width=4),
        name="Weights",
        showlegend=False
    )

    fig = tools.make_subplots(rows=3, cols=1, specs=[[{'rowspan': 2}],[None],[{}]])
    fig['layout'].update(height=600, width=800)
    fig['layout']['yaxis1'].update(title='$n_c$', range=[0,1], titlefont=dict(size=28), tickfont=dict(size=24))
    fig['layout']['yaxis2'].update(title='$w$', titlefont=dict(size=28), tickfont=dict(size=24))
    fig['layout']['xaxis2'].update(title='Percentage of slices', titlefont=dict(size=28), tickfont=dict(size=24))
    fig['layout']['xaxis1'].update(tickfont=dict(size=24))
    fig['layout']['font'].update(size=28)

    fig.append_trace(trace1,1,1)
    fig.append_trace(trace2,1,1)
    fig.append_trace(trace3,3,1)

    py.image.save_as(fig, 'paper2/'+filename+'.png')
    # plot_url = py.plot(fig, filename='basic-line')
#     print plot_url

def plot_graph_results(inc, inv):
    network = list(reversed(['intl', 'us', 'celegans', 'Chesapeake', 'ChesLower', 'ChesMiddle', 'ChesUpper', 'CrystalC', 'CrystalD', 'Everglades', 'Florida', 'Maspalomas', 'Michigan', 'Mondego', 'Narragan', 'Rhode', 'StMarks', 'one_mode_char', 'one_mode_message', 'days1-12', 'days1-4', 'days5-8', 'days9-12']))

    trace0 = go.Scatter(
        x=inc,
        y=network,
        mode='markers',
        name='Control distance for *-in_out-inc',
        marker=dict(
            color='#45A9DE',
            line=dict(
                color='#45A9DE',
                width=1,
            ),
            symbol='circle',
            size=16,
        )
    )
    trace1 = go.Scatter(
        x=inv,
        y=network,
        mode='markers',
        name='Control distance for *-in_out-inv',
        marker=dict(
            color='#BFCE80',
            line=dict(
                color='#BFCE80',
                width=1,
            ),
            symbol='circle',
            size=16,
        )
    )
    data = [trace0, trace1]
    layout = go.Layout(
        #title="Control distance for reassigned weights networks",
        font=dict(size=16),
        xaxis=dict(
            showgrid=False,
            showline=True,
            linecolor='rgb(102, 102, 102)',
            titlefont=dict(
                color='rgb(204, 204, 204)'
            ),
            tickfont=dict(
                color='rgb(102, 102, 102)',
            ),
            autotick=False,
            dtick=10,
            ticks='outside',
            tickcolor='rgb(102, 102, 102)',
        ),
        yaxis=dict(
            titlefont=dict(size=8)
        ),
        margin=dict(
            l=140,
            r=40,
            b=50,
            t=80
        ),
        legend=dict(
            font=dict(
                size=14,
            ),
            yanchor='middle',
            xanchor='right',
        ),
        width=800,
        height=600,
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        hovermode='closest',
    )
    fig = go.Figure(data=data, layout=layout)
    plot_url = py.plot(fig, filename='basic-line')
    print plot_url

def scatter_ly(r_er, er, r_ba, ba, r_ce, ce, r_us, us):
    trace0 = go.Scatter(
        x = r_us,
        y = us,
        mode='markers',
        name="us",
        marker=dict(
            color='rgba(254,209,82,0.9)',
            line=dict(
                color='rgba(254,209,82,1.0)',
                width=1,
            ),
            symbol='circle',
            size=16,
        )
    )

    trace1 = go.Scatter(
        x = r_er,
        y = er,
        mode='markers',
        name = "Erdos-Renyi",
        marker=dict(
            color='rgba(69,169,222,0.9)',
            line=dict(
                color='rgba(69,169,222,1.0)',
                width=1,
            ),
            symbol='circle',
            size=16,
        ),
    )

    trace2 = go.Scatter(
        x = r_ba,
        y = ba,
        mode='markers',
        name = "Barabasi-Albert",
        marker=dict(
            color='rgba(191,206,128,0.9)',
            line=dict(
                color='rgba(191,206,128,1.0)',
                width=1,
            ),
            symbol='circle',
            size=16,
        ),
    )

    trace3 = go.Scatter(
        x = r_ce,
        y = ce,
        mode='markers',
        name="celegans",
        marker=dict(
            color='rgba(240,84,89,0.9)',
            line=dict(
                color='rgba(240,84,89,1.0)',
                width=1,
            ),
            symbol='circle',
            size=16,
        )
    )

    data = [trace1, trace2, trace3, trace0]
    layout = go.Layout(
        font=dict(size=16),
        xaxis=dict(
            title="$\\rho$",
            showgrid=False,
            showline=True,
            linecolor='rgb(102, 102, 102)',
            titlefont=dict(
                size=28
            ),
            tickfont=dict(
            ),
            autotick=False,
            dtick=10,
            ticks='outside',
        ),
        yaxis=dict(
            title="Control distance",
            titlefont=dict(size=28)
        ),
        margin=dict(
            l=140,
            r=40,
            b=50,
            t=80
        ),
        legend=dict(
            font=dict(
                size=14,
            ),
            yanchor='middle',
            xanchor='middle',
        ),
        width=800,
        height=600,
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        hovermode='closest',
    )

    fig = go.Figure(data=data, layout=layout)
    plot_url = py.plot(fig, filename='basic-line')
    print plot_url
