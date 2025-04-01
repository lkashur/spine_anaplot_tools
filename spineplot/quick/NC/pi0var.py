# imports
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('pi0ana.mplstyle')
from scipy.optimize import curve_fit
from scipy import stats
import toml
import argparse

######################################################################
### Helper Functions
######################################################################
def get_stat_cov(bins):
    stat_cov = np.diag(bins)
    stat_err = np.sqrt(np.diagonal(stat_cov))
    return stat_cov, stat_err

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

def plot_histogram(cfg, plot_cfg, var, df_mc=None, pot_mc=0, df_onbeam=None, pot_onbeam=0, livetime_onbeam=0, df_offbeam=None, livetime_offbeam=0):
    ################
    ### Monte Carlo
    ################
    if plot_cfg['include_mc'] == True:
        contents_mc = []
        centers_mc = []
        labels_mc = []
        width_mc = []
        counts_mc = []
        frac_mc = []
        bins_mc, edges_mc = np.histogram(df_mc[var], bins=int(cfg['variables'][var]['bins'][0]), range=cfg['variables'][var]['bins'][1:])
        for i,m in enumerate(cfg['analysis_mc']['category_assignment']):
            mask = np.isin(df_mc[cfg['analysis_mc']['category_branch']], m)
            b, e = np.histogram(df_mc[var][mask], bins=int(cfg['variables'][var]['bins'][0]), range=cfg['variables'][var]['bins'][1:])
            contents_mc.append(b)
            centers_mc.append((e[1:] + e[:-1]) / 2.0)
            width_mc.append(np.diff(e))
            labels_mc.append(cfg['analysis_mc']['category_labels'][i])
            counts_mc.append(len(df_mc[var][mask]))
            frac_mc.append( len(df_mc[var][mask]) / len(df_mc[var]) )
        contents_mc = contents_mc[::-1]
        stat_cov_mc, stat_err_mc = get_stat_cov(bins_mc)
        err_mc = stat_err_mc
        if plot_cfg['normalization'] == 'data':
            norm_mc = np.outer(sum(contents_mc), sum(contents_mc))
            frac_stat_cov_mc = np.divide(stat_cov_mc, norm_mc, where=norm_mc!=0)
            bins_mc = (bins_mc / pot_mc) * pot_onbeam
            contents_mc = [(c / pot_mc) * pot_onbeam for c in contents_mc]
            counts_mc = [(c / pot_mc) * pot_onbeam for c in counts_mc]
            scaled_stat_cov_mc = np.outer(sum(contents_mc), sum(contents_mc)) * frac_stat_cov_mc            
            stat_err_mc = np.sqrt(np.diagonal(scaled_stat_cov_mc))
            err_mc = stat_err_mc

            #################
            ### Off-beam Data
            #################
            if plot_cfg['include_offbeam'] == True:
                bins_offbeam, edges_offbeam = np.histogram(df_offbeam[var], bins=int(cfg['variables'][var]['bins'][0]), range=cfg['variables'][var]['bins'][1:])
                centers_offbeam = (edges_offbeam[1:] + edges_offbeam[:-1]) / 2.0
                width_offbeam = np.diff(edges_offbeam)
                stat_cov_offbeam, stat_err_offbeam = get_stat_cov(bins_offbeam)
                norm_offbeam = np.outer(bins_offbeam, bins_offbeam)
                frac_stat_cov_offbeam = np.divide(stat_cov_offbeam, norm_offbeam, where=norm_offbeam!=0)
                bins_offbeam = (bins_offbeam / livetime_offbeam) * livetime_onbeam
                #print(sum(bins_offbeam))
                #print(np.sum(bins_offbeam))
                #print(livetime_onbeam / livetime_offbeam)
                
                ##################################################
                ### Add scaled off-beam content to cosmic category
                ##################################################
                contents_mc[0] = contents_mc[0] + bins_offbeam # 0
                counts_mc[-1] = counts_mc[-1] + sum(bins_offbeam)
                frac_mc = [c/sum(counts_mc) for c in counts_mc]
                
                #################################################
                ### Combine MC and off-beam error (stat. for now)
                #################################################
                scaled_stat_cov_offbeam = np.outer(bins_offbeam, bins_offbeam) * frac_stat_cov_offbeam
                scaled_stat_cov_mc_offbeam = scaled_stat_cov_mc + scaled_stat_cov_offbeam
                stat_err_mc = np.sqrt(np.diagonal(scaled_stat_cov_mc_offbeam))
                err_mc = stat_err_mc

        # back to MC...
        if plot_cfg['normalization'] == 'absolute':
            norm_mc = np.outer(sum(contents_mc), sum(contents_mc))
            frac_stat_cov_mc = np.divide(stat_cov_mc, norm_mc, where=norm_mc!=0)
            scale_mc = 1/np.sum([np.sum(i) for i in contents_mc])
            bins_mc = scale_mc * bins_mc
            contents_mc = [scale_mc * c for c in contents_mc]
            scaled_stat_cov_mc = np.outer(sum(contents_mc), sum(contents_mc)) * frac_stat_cov_mc
            stat_err_mc = np.sqrt(np.diagonal(scaled_stat_cov_mc))
            err_mc = stat_err_mc

    ################
    ### On-beam data
    ################
    if plot_cfg['include_onbeam'] == True:
            bins_onbeam, edges_onbeam = np.histogram(df_onbeam[var], bins=int(cfg['variables'][var]['bins'][0]), range=cfg['variables'][var]['bins'][1:])
            centers_onbeam = (edges_onbeam[1:] + edges_onbeam[:-1]) / 2.0
            width_onbeam = np.diff(edges_onbeam)
            label_onbeam = f'Data ({np.sum(bins_onbeam):.0f})'
            stat_cov_onbeam, stat_err_onbeam = get_stat_cov(bins_onbeam)
            if plot_cfg['normalization'] == 'absolute':
                frac_stat_cov_onbeam = stat_cov_onbeam / np.outer(bins_onbeam, bins_onbeam)
                scale_onbeam = 1/np.sum(bins_onbeam)
                bins_onbeam = scale_onbeam * bins_onbeam
                scaled_stat_cov_onbeam = np.outer(bins_onbeam, bins_onbeam) * frac_stat_cov_onbeam
                stat_err_onbeam  = np.sqrt(np.diag(scaled_stat_cov_onbeam))

    # Plot
    fig, ax = plt.subplots(figsize=(10,7))
    if plot_cfg['include_mc'] == True:
        colors = cfg['analysis_mc']['category_colors'][::-1]
        ax.hist(centers_mc, weights=contents_mc, bins=int(cfg['variables'][var]['bins'][0]), range=cfg['variables'][var]['bins'][1:], label=labels_mc, color=colors, histtype='barstacked')
        if plot_cfg['include_mc_err'] == True:
            ax.bar(x=centers_mc[0], height=2*stat_err_mc, bottom=np.sum(contents_mc, axis=0) - err_mc, width=width_mc[0], color='grey', hatch='xx', ec='darkgrey', alpha=0.4, lw=0.25, label='MC Uncertainty (Stat.)')
        

    if plot_cfg['include_onbeam'] == True:
        ax.errorbar(centers_onbeam, bins_onbeam, xerr=width_onbeam/2, yerr=stat_err_onbeam, linestyle='none', ecolor='black', label=label_onbeam)

    h, l = ax.get_legend_handles_labels()
    if plot_cfg['include_mc'] == True and plot_cfg['include_mc_err'] == True and plot_cfg['include_onbeam'] == False:
        h_mc, l_mc = h[:-1], l[:-1]
        l_mc  = [f'{l} ({counts_mc[li]:.0f}, {frac_mc[li]:.02%} )' for li,l in enumerate(l_mc)]
        #l_mc = [f'{l} ({np.sum(contents_mc[::-1][li]):.0f}, {np.sum(contents_mc[::-1][li]) / np.sum(contents_mc):.02%})'for li, l in enumerate(l_mc)]
        h_mc = h_mc[::-1]
        h_syst, l_syst = [h[-1]], [l[-1]]
        h, l = h_mc + h_syst, l_mc + l_syst
    elif plot_cfg['include_mc'] == True and plot_cfg['include_mc_err'] == True and plot_cfg['include_onbeam'] == True:
        h_mc, l_mc = h[:-2], l[:-2]
        l_mc  = [f'{l} ({counts_mc[li]:.0f}, {frac_mc[li]:.02%} )' for li,l in enumerate(l_mc)]
        #l_mc = [f'{l} ({np.sum(contents_mc[::-1][li]):.0f}, {np.sum(contents_mc[::-1][li]) / np.sum(contents_mc):.02%})'for li, l in enumerate(l_mc)]
        h_mc = h_mc[::-1]
        h_syst, l_syst = [h[-2]], [l[-2]]
        h_onbeam, l_onbeam = [h[-1]], [l[-1]]
        h, l = h_mc + h_syst + h_onbeam, l_mc + l_syst + l_onbeam

    ax.legend(h,l, ncols=2)
    if plot_cfg['normalization'] == 'data':
        pot_label = pot_onbeam
    else:
        pot_label = pot_mc
    add_plot_labels(ax, pot_label, adj_y=0.025, title=cfg['variables'][var]['title'])
    ax.set_xlim(([cfg['variables'][var]['bins'][1], cfg['variables'][var]['bins'][2]]))
    ax.set_xlabel(cfg['variables'][var]['xlabel'])
    ax.set_ylabel(cfg['variables'][var]['ylabel'])
    ax.set_ylim(bottom=0)
    plt.savefig(f'{var}.png')
        
######################################################################
### Main Analysis
######################################################################
def main(args):

    # Load pi0 config
    pi0_config = toml.load(args.config)

    # Load input file
    rf = uproot.open(args.in_file)
    
    # Load trees
    sel_nu_mc_tree = rf['events/mc/SelectedNu_NCCuts']
    sel_nu_mc_df = sel_nu_mc_tree.arrays(library='pd')
    sel_cos_mc_tree = rf['events/mc/SelectedCos_NCCuts']
    sel_cos_mc_df = sel_cos_mc_tree.arrays(library='pd')
    sel_mc_df = pd.concat([sel_nu_mc_df, sel_cos_mc_df])
    pot_mc = rf['events/mc/POT'].to_numpy()[0][0]
    livetime_mc = rf['events/mc/Livetime'].to_numpy()[0][0]

    #sel_offbeam_tree = rf['events/onbeam/SelectedCos_NCCuts']
    #sel_offbeam_df = sel_offbeam_tree.arrays(library='pd')
    #pot_offbeam = rf['events/offbeam/POT'].to_numpy()[0][0]
    #livetime_offbeam = rf['events/offbeam/Livetime'].to_numpy()[0][0]

    #sel_onbeam_tree = rf['events/onbeam/SelectedNu_NCCuts']
    #sel_onbeam_df = sel_onbeam_tree.arrays(library='pd')
    #pot_onbeam = rf['events/onbeam/POT'].to_numpy()[0][0]
    #livetime_onbeam = rf['events/onbeam/Livetime'].to_numpy()[0][0]

    # Plotting
    
    plot_config = {'include_mc' : True,
                   'include_mc_err': True,
                   'include_onbeam' : False,
                   'include_offbeam' : False,
                   'normalization' : None}
    
    for var in pi0_config['variables']:
        #plot_histogram(pi0_config, plot_config, var, sel_mc_df, pot_mc, sel_onbeam_df, pot_onbeam, livetime_onbeam, sel_offbeam_df, livetime_offbeam)
        #plot_histogram(pi0_config, plot_config, var, sel_mc_df, pot_mc, sel_onbeam_df, pot_onbeam)
        plot_histogram(pi0_config, plot_config, var, sel_mc_df, pot_mc)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--in_file', required=True)
    args = parser.parse_args()
    main(args)
