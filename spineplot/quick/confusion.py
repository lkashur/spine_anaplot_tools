# imports
import uproot
import ast
from ROOT import TFile, TEfficiency, TH1D, TGraphAsymmErrors, RDataFrame, TCanvas
import numpy as np
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
plt.style.use('pi0ana.mplstyle')
import pandas as pd
import toml
import argparse
from array import array
from evaluation import *

def main(args):
    
    # Load input
    rf = uproot.open(args.input)
    pot = rf['events/mc/POT'].to_numpy()[0][0]
    eff_pid_confusion_tree = rf['events/mc/EffPIDConfusion_PhaseCuts']
    eff_pid_confusion_df = eff_pid_confusion_tree.arrays(library="pd")
    #pur_pid_confusion_tree = rf['events/mc/PurPIDConfusion_PhaseCuts'] 

    # Plot
    plot_pid_confusion(eff_pid_confusion_df, pot, pur_or_eff='eff')

def add_plot_labels(ax, pot, adj_y=0.025, title=str()):
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.025, #0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[0] + 0.05*(xrange[1] - xrange[0]) #0.025
    prelim_label = r'$\bf{ICARUS}$ Work-in-Progress'
    ax.text(x=usex, y=usey, s=prelim_label, fontsize=14, color='#d67a11')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = (xrange[1] + xrange[0]) / 2
    ax.text(x=usex, y=usey, s=f'{title}', fontsize=14, fontweight='bold', color='black')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[1] - 0.05*(xrange[1] - xrange[0])
    mag = int(np.floor(np.log10(pot)))
    usepot = pot/10**mag
    s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')

def plot_pid_confusion(pid_confusion_df, pot, pur_or_eff='eff'):
    # Start from true particles for efficiency; reco particles for purity                                                      
    if(pur_or_eff == 'eff'):
        _from = 'true'
        _to = 'reco'
    if(pur_or_eff == 'pur'):
        _from = 'reco'
        _to = 'true'

    pid_confusion_df.dropna(inplace=True)

    # Get particles
    photons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 0]
    photon_pids_from = [0] * len(photons)
    photon_pids_to = photons[f'{_to}_pid'].to_list()

    electrons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 1]
    electron_pids_from = [1] * len(electrons)
    electron_pids_to = electrons[f'{_to}_pid'].to_list()

    muons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 2]
    muon_pids_from = [2] * len(muons)
    muon_pids_to = muons[f'{_to}_pid'].to_list()

    pions = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 3]
    pion_pids_from = [3] * len(pions)
    pion_pids_to = pions[f'{_to}_pid'].to_list()

    protons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 4]
    proton_pids_from = [4] * len(protons)
    proton_pids_to = protons[f'{_to}_pid'].to_list()

    kaons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 5]
    kaon_pids_from = [5] * len(kaons)
    kaon_pids_to = kaons[f'{_to}_pid'].to_list()

    # One flat list for all particles                                             
    pids_from = photon_pids_from + electron_pids_from + muon_pids_from + pion_pids_from + proton_pids_from + kaon_pids_from
    pids_to = photon_pids_to + electron_pids_to + muon_pids_to + pion_pids_to + proton_pids_to + kaon_pids_to

    # Set labels
    label_keys = np.unique(pids_from)
    label_keys = np.append(label_keys, 3) # show pions                                                                                      
    label_keys = sorted(label_keys)
    label_values = []
    if 0 in label_keys: label_values.append('Photon')
    if 1 in label_keys: label_values.append('Electron')
    if 2 in label_keys: label_values.append('Muon')
    #if 3 in label_keys: label_values.append('Pion')
    label_values.append('Pion') # show pions
    if 4 in label_keys: label_values.append('Proton')
    if 5 in label_keys: label_values.append('Kaon')

    # Plot                    
    cm = confusion_matrix(np.array(pids_from), np.array(pids_to), normalize='true', labels=label_keys)
    cm_counts = confusion_matrix(np.array(pids_from), np.array(pids_to), labels=label_keys)
    fig, ax = plt.subplots(figsize=(8, 8))

    im = heatmap(cm.T * 100,
                 label_values,
                 label_values,
                 ax=ax, cmap="Blues")
    annotate_heatmap(im, unc=cm_counts.T, valfmt="{:.2f}% \n({})", fontsize=14)
    plt.setp(ax.get_yticklabels(), rotation=90, ha="center",rotation_mode="anchor")
    ax.tick_params(axis='both', which='major', pad=15, bottom=False, top=False, left=False, right=False, labelsize=14)
    ax.tick_params(axis="x", pad=1)
    if(pur_or_eff == 'eff'): ax.set_xlabel('True Signal Particles', fontsize=15, labelpad=10)
    if(pur_or_eff == 'pur'): ax.set_xlabel('Selected Reco Particles', fontsize=15, labelpad=10)
    ax.xaxis.set_label_position('bottom')
    if(pur_or_eff == 'eff'): ax.set_ylabel('Matched Reco Particles', fontsize=15)
    if(pur_or_eff == 'pur'): ax.set_ylabel('Matched True Particles', fontsize=15)
    add_plot_labels(ax,pot, adj_y=0.045, title=str())
    plt.savefig(f'{pur_or_eff}_pid_confusion.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    args = parser.parse_args()
    main(args)
