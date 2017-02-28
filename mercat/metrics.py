#!/usr/bin/env python

import plotly.graph_objs as go
from plotly.offline import plot


#predict_isoelectric_point_ProMoST code from
#http://isoelectric.ovh.org/
#http://biologydirect.biomedcentral.com/articles/10.1186/s13062-016-0159-9
#IPC - Isoelectric Point Calculator
#Lukasz P. Kozlowski
#Biology Direct 2016
#DOI: 10.1186/s13062-016-0159-9


# N-teminus, middle, C-terminus
promost = {
    'K': [10.00, 9.80, 10.30],
    'R': [11.50, 12.50, 11.50],
    'H': [4.89, 6.08, 6.89],
    'D': [3.57, 4.07, 4.57],
    'E': [4.15, 4.45, 4.75],
    'C': [8.00, 8.28, 9.00],
    'Y': [9.34, 9.84, 10.34],
    'U': [5.20, 5.43, 5.60],  # ref (http://onlinelibrary.wiley.com/doi/10.1002/bip.21581/pdf)
}

promost_mid = {
    "G": [7.50, 3.70],
    "A": [7.58, 3.75],
    "S": [6.86, 3.61],
    "P": [8.36, 3.40],
    "V": [7.44, 3.69],
    "T": [7.02, 3.57],
    "C": [8.12, 3.10],
    "I": [7.48, 3.72],
    "L": [7.46, 3.73],
    "J": [7.46, 3.73],
    "N": [7.22, 3.64],
    "D": [7.70, 3.50],
    "Q": [6.73, 3.57],
    "K": [6.67, 3.40],
    "E": [7.19, 3.50],
    "M": [6.98, 3.68],
    "H": [7.18, 3.17],
    "F": [6.96, 3.98],
    "R": [6.76, 3.41],
    "Y": [6.83, 3.60],
    "W": [7.11, 3.78],
    "X": [7.26, 3.57],  # avg
    "Z": [6.96, 3.535],  # ("E"+"Q")/2
    'B': [7.46, 3.57],  # ("N"+"D")/2
    'U': [5.20, 5.60],
    'O': [7.00, 3.50],
}

import sys
def predict_isoelectric_point_ProMoST(seq):
    '''Calculate isoelectric point using ProMoST model'''
    NQ = 0.0
    pH = 6.51  # starting po pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability of finding the solution
    pHprev = 0.0
    pHnext = 14.0
    E = 0.01  # epsilon means precision [pI = pH +- E]
    temp = 0.01
    while 1:
        if seq[0] in promost:
            QN1 = -1.0 / (1.0 + pow(10, (promost[seq[0]][2] - pH)))
        else:
            QN1 = -1.0 / (1.0 + pow(10, (promost_mid[seq[0]][1] - pH)))
        # print
        if seq[-1] in promost:
            QP2 = 1.0 / (1.0 + pow(10, (pH - promost[seq[-1]][0])))
        elif seq[-1] in promost_mid:
            QP2 = 1.0 / (1.0 + pow(10, (pH - promost_mid[seq[-1]][0])))
        else:
            print(seq[-1] + " not found!")

            sys.exit(1)

        QN2 = -seq.count('D') / (1.0 + pow(10, (promost['D'][1] - pH)))
        QN3 = -seq.count('E') / (1.0 + pow(10, (promost['E'][1] - pH)))
        QN4 = -seq.count('C') / (1.0 + pow(10, (promost['C'][1] - pH)))
        QN5 = -seq.count('Y') / (1.0 + pow(10, (promost['Y'][1] - pH)))
        QP1 = seq.count('H') / (1.0 + pow(10, (pH - promost['H'][1])))
        QP3 = seq.count('K') / (1.0 + pow(10, (pH - promost['K'][1])))
        QP4 = seq.count('R') / (1.0 + pow(10, (pH - promost['R'][1])))

        NQ = QN1 + QN2 + QN3 + QN4 + QN5 + QP1 + QP2 + QP3 + QP4
        # %%%%%%%%%%%%%%%%%%%%%%%%%   BISECTION   %%%%%%%%%%%%%%%%%%%%%%%%
        if NQ < 0.0:  # we are out of range, thus the new pH value must be smaller
            temp = pH
            pH = pH - ((pH - pHprev) / 2.0)
            pHnext = temp
            # print "pH: ", pH, ", \tpHnext: ",pHnext
        else:
            temp = pH
            pH = pH + ((pHnext - pH) / 2.0)
            pHprev = temp
            # print "pH: ", pH, ",\tpHprev: ", pHprev

        if (pH - pHprev < E) and (pHnext - pH < E):  # terminal condition, finding pI with given precision
            return round(pH,2)


mass_aa = { #monoisotopic 	average #Kyte-Doolittle?
"A"	: 71.0788,
"B" : 114.6686,   # Asx   Aspartic acid or Asparagine
"C"	: 103.1388,
"D"	: 115.0886,
"E"	: 129.1155,
"F"	: 147.1766,
"G"	: 57.0519,
"H"	: 137.1411,
"I"	: 113.1594,
"K"	: 128.1741,
"L"	: 113.1594,
"M"	: 131.1926,
"N"	: 114.1038,
"O"	: 237.3018,
"P"	: 97.1167,
"Q"	: 128.1307,
"R"	: 156.1875,
"S"	: 87.0782,
"T"	: 101.1051,
"U"	: 150.0388,
"V"	: 99.1326,
"W"	: 186.2132,
"X" : 111.1138, # Xaa   Any amino acid
"Y"	: 163.176,
"Z" : 128.7531    #Glx   Glutamine or Glutamic acid
}

# Kyte-Doolittle Hydropathy Scores
hydro_scores = {
"A" : 1.8,
"R" : -4.5,
"N" : -3.5,
"D" : -3.5,
"C" : 2.5,
"Q" : -3.5,
"E" : -3.5,
"G" : -0.4,
"H" : -3.2,
"I" : 4.5,
"L" : 3.8,
"K" : -3.9,
"M" : 1.9,
"F" : 2.8,
"P" : -1.6,
"S" : -0.8,
"T" : -0.7,
"W" : -0.9,
"Y" : -1.3,
"V" : 4.2

 }


def calculate_MW(seq): #calculate based on average mass
    mw = 0.0
    for c in seq:
        if c in mass_aa: mw += mass_aa[c]
    mw += 18.01524 #avg for water molecule
    return round(mw,2)


def calculate_hydro(seq):
    hydro = 0.0
    for c in seq:
        if c in hydro_scores: hydro += hydro_scores[c]
    return round(hydro,2)


from skbio.diversity import alpha as skbio_alpha
from skbio.diversity import beta as skbio_beta

def mercat_compute_alpha_beta_diversity(counts,bif):

    abm = dict()

    abm['shannon'] = skbio_alpha.shannon(counts)
    abm['simpson'] = skbio_alpha.simpson(counts)
    abm['simpson_e'] = skbio_alpha.simpson_e(counts)
    abm['goods_coverage'] = skbio_alpha.goods_coverage(counts)
    abm['fisher_alpha'] = skbio_alpha.fisher_alpha(counts)
    abm['dominance'] = skbio_alpha.dominance(counts)
    abm['chao1'] = skbio_alpha.chao1(counts)
    abm['chao1_ci'] = skbio_alpha.chao1_ci(counts)
    abm['ace'] = skbio_alpha.ace(counts)

    with open(bif + "_diversity_metrics.txt", 'w') as dmptr:
        for abmetric in abm:
            dmptr.write(abmetric + " = " + str(abm[abmetric]) + "\n")



def mercat_scatter_plots(bif,xlab,res_df,kmerstring):

    axis_title_font_size = 20
    axis_tick_label_size = 18
    legend_font_size = 14
    marker_size = 10

    save_img = 0

    trace1 = go.Scatter(
        x=res_df[xlab].values,
        y=res_df.index.values,
        mode='markers'
    )

    data = go.Data([trace1])
    layout = go.Layout(
        legend=dict(
            font=dict(
                # family='sans-serif',
                size=legend_font_size,
                color='black'
            ),
            # borderwidth=2
        ),
        autosize=True,
        height=632,
        # title='Baseline',
        width=1274,
        xaxis=dict(
            # autorange=True,
            fixedrange=False,
            title=xlab,
            # type='linear',
            # showgrid=False,
            ticks='inside',
            ticklen=8,
            tickwidth=2,
            tickcolor='#000',
            #tickvals=list(reversed(cores)),
            titlefont=dict(
                # family='Courier New, monospace',
                size=axis_title_font_size,
                color='black'
            ),
            tickfont=dict(
                # family='Old Standard TT, serif',
                size=axis_tick_label_size,
                color='black'
            ),
            showgrid=True,
            showline=True,
            mirror=True,
            zeroline=False,
            # gridcolor='black',
            # linecolor='black',
            gridwidth=2
        ),
        yaxis=dict(
            #type='log',
            autorange=True,
            fixedrange=False,
            ticks='inside',
            ticklen=8,
            tickwidth=2,
            tickcolor='#000',
            titlefont=dict(
                # family='Courier New, monospace',
                size=axis_title_font_size,
                color='black'
                # color='#7f7f7f'
            ),
            rangemode='normal',
            tickmode='linear',
            # tickwidth=4,
            tickfont=dict(
                # family='Old Standard TT, serif',
                size=axis_tick_label_size,
                color='black'
            ),
            showgrid=True,
            showline=True,
            mirror=True,
            zeroline=False,
            # gridcolor='black',
            # linecolor='black',
            gridwidth=2,
            title=kmerstring
        )
    )

    fig = go.Figure(data=data, layout=layout)
    plot(fig, filename=bif + "_"+xlab+".html", auto_open=False)



def mercat_stackedbar_plots(bif,xlab,res_df,kmerstring,total_freq_count):

    axis_title_font_size = 20
    axis_tick_label_size = 18
    legend_font_size = 14
    marker_size = 10

    #total_freq_count = res_df[xlab].sum()
    index_vals = res_df.index.values

    btraces = []
    topk10 = min(len(index_vals),10)
    for i in range(0,topk10):
        fr = res_df.loc[index_vals[i], 'Count'],
        trace1 = go.Bar(
            x=[bif],
            y=[(fr[0]*1.0/total_freq_count)],
            name = index_vals[i]
        )
        btraces.append(trace1)

    data = go.Data(btraces)
    layout = go.Layout(
        barmode='stack',
        legend=dict(
            font=dict(
                # family='sans-serif',
                size=legend_font_size,
                color='black'
            ),
            # borderwidth=2
        ),
        autosize=True,
        height=632,
        # title='Baseline',
        width=1274,
        xaxis=dict(
            # autorange=True,
            fixedrange=False,
            title="Samples",
            # type='linear',
            # showgrid=False,
            ticks='inside',
            ticklen=8,
            tickwidth=2,
            tickcolor='#000',
            #tickvals=list(reversed(cores)),
            titlefont=dict(
                # family='Courier New, monospace',
                size=axis_title_font_size,
                color='black'
            ),
            tickfont=dict(
                # family='Old Standard TT, serif',
                size=axis_tick_label_size,
                color='black'
            ),
            showgrid=True,
            showline=True,
            mirror=True,
            zeroline=False,
            # gridcolor='black',
            # linecolor='black',
            gridwidth=2
        ),
        yaxis=dict(
            #type='log',
            autorange=True,
            fixedrange=False,
            ticks='inside',
            ticklen=8,
            tickwidth=2,
            tickcolor='#000',
            titlefont=dict(
                # family='Courier New, monospace',
                size=axis_title_font_size,
                color='black'
                # color='#7f7f7f'
            ),
            rangemode='normal',
            tickmode='linear',
            # tickwidth=4,
            tickfont=dict(
                # family='Old Standard TT, serif',
                size=axis_tick_label_size,
                color='black'
            ),
            showgrid=True,
            showline=True,
            mirror=True,
            zeroline=False,
            # gridcolor='black',
            # linecolor='black',
            gridwidth=2,
            title=kmerstring
        )
    )

    fig = go.Figure(data=data, layout=layout)
    plot(fig, filename=bif + "_barchart_"+xlab+".html", auto_open=False)