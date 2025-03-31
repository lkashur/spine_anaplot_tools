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
    pot = rf['events/mc/POT'].to_numpy()[0][0]
    muon_tree = rf['events/mc/Muon']
    muon_df = muon_tree.arrays(library="pd")
    photon_tree = rf['events/mc/Photon']
    photon_df = photon_tree.arrays(library="pd")

    # Limit each tree to only cases where true particle is matched to reco particle that is primary, contained, and has same PID
    muon_df = muon_df[(muon_df.true_is_contained == True) & (muon_df.reco_is_primary == True) & (muon_df.reco_is_contained == True) & (muon_df.reco_is_muon == True)]
    photon_df = photon_df[(photon_df.reco_is_primary == True) & (photon_df.reco_is_contained == True) & (photon_df.reco_is_photon == True)]

    #####################################################################################
    #####################################################################################

    '''
    muon_df['muon_ke_res'] = (muon_df.reco_calo_ke - muon_df.reco_csda_ke) / muon_df.reco_csda_ke
    plt.hist(muon_df['muon_ke_res'], bins=30, range=[-0.5, 0.5])
    plt.xlabel('(Calo. KE - CSDA KE) / CSDA KE')
    plt.ylabel('Counts')
    plt.show()
    '''

    plot_muon_energy_bias(muon_df, pot, 'true_ke', 'reco_csda_ke')
    #plot_muon_energy_bias(muon_df, pot, reco_energy='reco_mcs_ke')
    #plot_shower_corr(photon_df, pot)
    #plot_shower_energy_bias(photon_df, pot)


def plot_muon_energy_bias(muon_df, pot, x='reco_csda_ke', y='reco_calo_ke'):
    fig, ax = plt.subplots(figsize=(10,7))
    h = ax.hist2d(muon_df[x], muon_df[y], bins=np.linspace(0,2000, 101), cmap='viridis', norm=colors.LogNorm())
    ax.plot([0, 2000], [0, 2000], color='red', linestyle='--', lw=1.5)
    ax.set_xlabel('True Muon Kinetic Energy [MeV]')
    ax.set_ylabel('Reco Muon CSDA Kinetic Energy [MeV]')
    cbar = fig.colorbar(h[3], ax=ax)
    cbar.set_label('Entries')
    add_plot_labels(ax, pot, adj_y=0.020, title=str())
    plt.savefig(f'muon_{x}_vs_{y}.png')

'''    
def plot_muon_energy_bias(muon_df, pot, reco_energy='reco_csda_ke'):
    fig, ax = plt.subplots(figsize=(10,7))
    h = ax.hist2d(muon_df.true_muon_ke, muon_df[reco_energy], bins=np.linspace(0,2000, 201), cmap='viridis', norm=colors.LogNorm())
    ax.plot([0, 2000], [0, 2000], color='red', linestyle='--', lw=1.5)
    ax.set_xlabel('True Muon Kinetic Energy [MeV]')
    ax.set_ylabel('Reconstructed Muon Kinetic Energy [MeV]')
    cbar = fig.colorbar(h[3], ax=ax)
    cbar.set_label('Entries')
    add_plot_labels(ax, pot, adj_y=0.020, title=str())
    plt.savefig(f'muon_{reco_energy}_vs_true_ke.png')
'''

def plot_shower_energy_bias(photon_df, pot):
    fig, ax = plt.subplots(figsize=(10,7))
    h = ax.hist2d(photon_df.true_photon_ke, photon_df.reco_calo_ke, bins=np.linspace(0,1000, 51), cmap='viridis', norm=colors.LogNorm())
    ax.plot([0, 1000], [0, 1000], color='red', linestyle='--', lw=1.5)
    ax.set_xlabel('True Photon Energy [MeV]')
    ax.set_ylabel('Reconstructed Photon Energy [MeV]')
    cbar = fig.colorbar(h[3], ax=ax)
    cbar.set_label('Entries')
    add_plot_labels(ax, pot, adj_y=0.020, title=str())
    plt.savefig(f'photon_reco_vs_true_energy.png')

def plot_shower_corr(photon_df, pot):
    photon_df['frac'] = photon_df.reco_calo_ke_pre_corr / photon_df.true_photon_ke
    photon_df['e_true_bin'] = pd.cut(photon_df.true_photon_ke, np.linspace(0,1000,21))

    e_mid = []
    frac_agg = []
    frac_low = []
    frac_high = []
    for i,(name,group) in enumerate(photon_df.groupby('e_true_bin')):

        if len(group) < 2: continue

        frac_temp = (group.frac.to_list(),)
        bootstrap_ci = bootstrap(frac_temp, np.median, n_resamples=1000, method='percentile', confidence_level=0.68)

        e_mid.append(group['e_true_bin'].iloc[0].mid)
        frac_agg.append(np.median(group.frac))
        frac_low.append(bootstrap_ci.confidence_interval[0])
        frac_high.append(bootstrap_ci.confidence_interval[1])

    frac_agg = np.array(frac_agg)
    frac_low = np.array(frac_low)
    frac_high = np.array(frac_high)
    fudge = np.mean(frac_agg)

    fig,ax = plt.subplots(figsize=(10,7))
    h = ax.hist2d(x=photon_df.true_photon_ke, y=photon_df.frac, bins=(20,80), range=((0,1000),(0,1.2)), cmap='Blues', norm=colors.LogNorm())
    ax.axhline(fudge, c='limegreen', label='E$_{reco}$/E$_{true}$ = ' + str(round(fudge, 2)), zorder=1)
    ax.errorbar(e_mid, frac_agg, yerr=((frac_agg - frac_low),(frac_high - frac_agg)), ls='None', capsize=3, ecolor='orange')
    ax.scatter(e_mid, frac_agg, c='orange', s=10)

    ax.set_xlabel('True Photon Energy [MeV]')
    ax.set_ylabel('Completeness')
    ax.set_xlim((0,800))
    ax.set_ylim((0,1))
    ax.legend(loc='lower right', fontsize=14)
    cbar = fig.colorbar(h[3], ax=ax)
    cbar.set_label('Entries')
    add_plot_labels(ax, pot, adj_y=0.020, title=str())
    plt.savefig('shower_corr.png')

def add_plot_labels(ax, pot, adj_y=0.025, title=str()):
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.025, #0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[0] + 0.025*(xrange[1] - xrange[0]) #0.025
    prelim_label = r'$\bf{ICARUS}$ Simulation Work-in-Progress'
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
