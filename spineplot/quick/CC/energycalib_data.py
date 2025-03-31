# imports
import uproot
from ROOT import TFile, TEfficiency, TH1D, TGraphAsymmErrors, RDataFrame, TCanvas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
plt.style.use('pi0ana.mplstyle')
import pandas as pd
import toml
from scipy.stats import bootstrap
import argparse
from array import array

def main(args):

    # Load ROOT file into pandas dataframe
    rf = uproot.open(args.in_file)
    sel_onbeam_tree = rf['events/onbeam/SelectedNu_TradCuts']
    sel_onbeam_df = sel_onbeam_tree.arrays(library='pd')
    pot_onbeam = rf['events/onbeam/POT'].to_numpy()[0][0]
    

    #####################################################################################
    #####################################################################################

    
    sel_onbeam_df['muon_ke_res'] = (sel_onbeam_df.muon_calo_ke - sel_onbeam_df.muon_csda_ke) / sel_onbeam_df.muon_csda_ke
    plt.hist(sel_onbeam_df['muon_ke_res'], bins=30, range=[-0.5, 0.5])
    plt.xlabel('(Calo. KE - CSDA KE) / CSDA KE')
    plt.ylabel('Counts')
    plt.show()
    #plot_muon_energy_bias(sel_onbeam_df, pot_onbeam)
    
def plot_muon_energy_bias(sel_onbeam_df, pot):
    fig, ax = plt.subplots(figsize=(10,7))
    h = ax.hist2d(sel_onbeam_df.muon_csda_ke, sel_onbeam_df.muon_calo_ke, bins=np.linspace(0,1000, 101), cmap='viridis', cmin=1)
    ax.plot([0, 1000], [0, 1000], color='red', linestyle='--', lw=1.5)
    ax.set_xlabel('CSDA Muon Kinetic Energy [MeV]')
    ax.set_ylabel('Calo. Muon Kinetic Energy [MeV]')
    cbar = fig.colorbar(h[3], ax=ax)
    cbar.set_label('Entries')
    add_plot_labels(ax, pot, adj_y=0.020, title=str())
    plt.savefig(f'muon_csda_vs_calo_ke.png')

def add_plot_labels(ax, pot, adj_y=0.025, title=str()):
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.025, #0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[0] + 0.025*(xrange[1] - xrange[0]) #0.025
    prelim_label = r'$\bf{ICARUS}$ Work-in-Progress'
    ax.text(x=usex, y=usey, s=prelim_label, fontsize=14, color='#d67a11')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = (xrange[1] - xrange[0]) / 2
    ax.text(x=usex, y=usey, s=f'{title}', fontsize=14, color='black')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[1] - 0.025*(xrange[1] - xrange[0])
    mag = int(np.floor(np.log10(pot)))
    usepot = pot/10**mag
    s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', required=True)
    args = parser.parse_args()
    main(args)
