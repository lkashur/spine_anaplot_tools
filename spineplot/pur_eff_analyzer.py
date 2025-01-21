# imports
import uproot
from ROOT import TFile, TEfficiency, TH1D, TGraphAsymmErrors, RDataFrame, TCanvas
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('configurations/pi02024/styles/pi02024.mplstyle')
import pandas as pd
import toml
import argparse
from array import array

def main(args):
    
    # Load variable configuration
    var_cfg = toml.load(args.var_cfg)

    # Load ROOT file into pandas dataframe
    rf = uproot.open(args.in_file)
    sel_tree = rf['events/mc/Selected_TradCuts']
    sel_df = sel_tree.arrays(library="pd")
    sig_tree = rf['events/mc/Signal_TradCuts']
    sig_df = sig_tree.arrays(library="pd")
    sig_df = sig_df[sig_df.category == 0] # only fiducialized signal events

    # Purity
    purity = calc_flat_purity(sel_df)
    print(f'Purity: {round(purity,2)}%')

    # Efficiency
    efficiency = calc_flat_efficiency(sig_df)
    print(f'Efficiency: {round(efficiency,2)}%')

    # Calculate efficiency as function of...
    #plot_diff_pur_eff(sig_df, 'muon_momentum_mag', var_cfg, [0,3200], 40, 'eff') # muon momentum
    #plot_diff_pur_eff(sig_df, 'pi0_momentum_mag', var_cfg, [0,1000], 30, 'eff') # pi0 momentum

def calc_flat_purity(sel_df):
    total_selected_events = len(sel_df)
    matched_selected_events = len(sel_df[sel_df.category == 0])
    return 100 * matched_selected_events / total_selected_events

def calc_flat_efficiency(sig_df):
    total_signal_events = len(sig_df)
    selected_signal_events = len(sig_df[sig_df['all_cut'] == 1])
    return 100 * selected_signal_events / total_signal_events

def plot_diff_pur_eff(df, var, var_cfg, var_range, nbins, metric):
    bcs = []
    bxerr0s = []
    bxerr1s = []
    byerr0s = []
    byerr1s = []
    vals = []

    # filtering/cuts                                                                                      
    if var_range[0] == var_range[-1]:
        df = df
    else:
        df = df[(df[var] > var_range[0]) & (df[var] < var_range[1])]

    #df['var_q'] = pd.qcut(df[var], q=nbins)
    df['var_q'] = pd.cut(df[var], np.linspace(var_range[0], var_range[1], nbins+1))

    hedges = sorted([i.left for i in df.var_q.unique().tolist()])
    lastedge = max([i.right for i in df.var_q.unique().tolist()])
    hedges.append(lastedge)
    hedges = np.array(hedges)
    hpass = TH1D('hpass', '', nbins, hedges)
    htotal = TH1D('htotal', '', nbins, hedges)
    counts = []

    for i,(name,group) in enumerate(df.groupby('var_q')):
        bxerr0s.append(name.left)
        bcs.append(name.mid)
        bxerr1s.append(name.right)
        if metric == 'pur':
            vals.append(calc_flat_purity(group) / 100)
            hpass.SetBinContent(i+1, len(group[group['category'] == 0]))
        elif metric == 'eff':
            vals.append(calc_flat_efficiency(group) / 100)
            hpass.SetBinContent(i+1, len(group[group['all_cut'] == 1]))
        htotal.SetBinContent(i+1, len(group))
        counts.append(len(group))

    gr = TGraphAsymmErrors()
    gr.Divide(hpass, htotal, 'cl=0.683 b(1,1) mode')                                             
    for i in range(nbins):
        byerr0s.append(gr.GetErrorYlow(i))
        byerr1s.append(gr.GetErrorYhigh(i))

    fig, ax1 = plt.subplots(figsize=(10,6))
    ax1.bar(x=hedges[:-1], height=counts, width=np.diff(hedges), align='edge', fc='C1', ec='black')
    ax1.set_xlabel(var_cfg['variables'][var]['xlabel'])
    ax1.set_ylabel('Entries')
    ax2 = ax1.twinx()
    ax2.set_ylim([0,1])
    ax2.errorbar(bcs, vals, xerr=[np.array(bcs) - np.array(bxerr0s), np.array(bxerr1s) - np.array(bcs)], yerr=[np.array(byerr0s), np.array(byerr1s)], fmt='o', capsize=2)
    if metric == 'pur':
        ylabel = 'Purity'
    elif metric == 'eff':
        ylabel = 'Efficiency'
    ax2.set_ylabel(ylabel)
    ax2.spines['right'].set_color('C0')
    ax2.yaxis.label.set_color('C0')
    ax2.tick_params(axis='y', colors='C0')
    plt.title('Signal $1\mu 0\pi^{±} 1\pi^{0}$')
    plt.savefig(metric + '_vs_' + var + '.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', required=True)
    parser.add_argument('--var_cfg', required=True)
    args = parser.parse_args()
    main(args)
