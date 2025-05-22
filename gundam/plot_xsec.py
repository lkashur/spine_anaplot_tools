# imports
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('pi0ana.mplstyle')
import argparse

def load_hist1d(filename, key):
    bins,edges = filename[key].to_numpy()
    centers = (edges[1:] + edges[:-1]) / 2.0
    errs = filename[key].errors()
    binwidth = np.diff(edges)
    return bins, edges, centers, errs, binwidth

def add_plot_labels(ax, pot, adj_y=0.025, title=str()):
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.025, #0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[0] + 0.05*(xrange[1] - xrange[0]) #0.025 
    prelim_label = r'$\bf{ICARUS}$ Work-in-Progress'
    #ax.text(x=usex, y=usey, s=prelim_label, fontsize=14, color='#d67a11')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
 
    xrange = ax.get_xlim()
    usex = (xrange[1] + xrange[0]) / 2
    ax.text(x=usex, y=usey, s=f'{title}', fontsize=14, fontweight='bold', color='black')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix

    #xrange = ax.get_xlim()
    usex = xrange[1] - 0.1*(xrange[1] - xrange[0])
    mag = int(np.floor(np.log10(pot)))
    usepot = pot/10**mag
    s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')

def main(args):
    
    # input
    rf = uproot.open(args.input)
    pot = rf['events/mc/POT'].to_numpy()[0][0]

    # config
    xlabel = r'$\mathrm{cos}(\mathrm{\theta_{\mu}})$'
    ylabel = r'$\frac{\mathrm{d\sigma}}{\mathrm{dcos(\theta_{\mu})}} [\frac{\mathrm{cm^{2}}}{\mathrm{nucleon}}]$'

    # load asimov
    rf_asimov = uproot.open(args.asimov)
    bins_asimov, edges_asimov, centers_asimov, errors_asimov, binwidth_asimov = load_hist1d(rf_asimov, 'calcXsec/histograms/TrueSignal_TrueMuonCos_TH1D')    
    # replace bins
    bins_asimov_temp, edges_asimov, centers_asimov, errors_asimov_temp, binwidth_asimov = load_hist1d(rf_asimov, 'calcXsec/plots/histograms/TrueSignal_TrueMuonCos/true_muon_beam_costheta/Data_TH1D')
    
    # load data
    rf_data = uproot.open(args.data)
    bins_data, edges_data, centers_data, errors_data, binwidth_data = load_hist1d(rf_asimov, 'calcXsec/histograms/TrueSignal_TrueMuonCos_TH1D')
    #bins_data = bins_data * 1.5

    # plot
    fig, ax = plt.subplots(figsize=(8,8))
    ax.errorbar(centers_asimov, bins_data, yerr=errors_data, fmt='o', color='black', ecolor='black', capsize=3, label='Run 2 Data', zorder=100)
    ax.plot(edges_asimov, np.insert(bins_asimov,-1, bins_asimov[-1]),c='red', linestyle='--', drawstyle='steps-post', label='GENIE AR23_20i')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim([-1, 1])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])
    add_plot_labels(ax,pot, adj_y=0.040, title=str())
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--asimov', required=True)
    parser.add_argument('--data', required=True)
    args = parser.parse_args()
    main(args)
