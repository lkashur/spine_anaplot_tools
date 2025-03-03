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

def crystalball(x, a, b, m, l, s, c, d):
    return a*s*stats.crystalball.pdf(x, b, m, l, s)+ c*x + d

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
    ax.text(x=usex, y=usey, s=f'{title}', fontsize=14, fontweight='bold', color='black')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[1] - 0.025*(xrange[1] - xrange[0])
    mag = int(np.floor(np.log10(pot)))
    usepot = pot/10**mag
    s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')

def plot_histogram(cfg, plot_cfg, var, df_mc=None, pot_mc=0, df_onbeam=None, pot_onbeam=0):
    ################
    ### Monte Carlo
    ################
    if plot_cfg['include_mc'] == True:
        contents_mc = []
        centers_mc = []
        labels_mc = []
        width_mc = []
        bins_mc, edges_mc = np.histogram(df_mc[var], bins=int(cfg['variables'][var]['bins'][0]), range=cfg['variables'][var]['bins'][1:])
        for i,m in enumerate(cfg['analysis_mc']['category_assignment']):
            mask = np.isin(df_mc[cfg['analysis_mc']['category_branch']], m)
            b, e = np.histogram(df_mc[var][mask], bins=int(cfg['variables'][var]['bins'][0]), range=cfg['variables'][var]['bins'][1:])
            contents_mc.append(b)
            centers_mc.append((e[1:] + e[:-1]) / 2.0)
            width_mc.append(np.diff(e))
            labels_mc.append(cfg['analysis_mc']['category_labels'][i])
        contents_mc = contents_mc[::-1]
        stat_cov_mc, stat_err_mc = get_stat_cov(bins_mc)
        err_mc = stat_err_mc
        if plot_cfg['normalization'] == 'data':
            norm_mc = np.outer(sum(contents_mc), sum(contents_mc))
            frac_stat_cov_mc = np.divide(stat_cov_mc, norm_mc, where=norm_mc!=0)
            bins_mc = (bins_mc / pot_mc) * pot_onbeam
            contents_mc = [(c / pot_mc) * pot_onbeam for c in contents_mc]
            scaled_stat_cov_mc = np.outer(sum(contents_mc), sum(contents_mc)) * frac_stat_cov_mc
            stat_err_mc = np.sqrt(np.diagonal(scaled_stat_cov_mc))
            err_mc = stat_err_mc
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
            ax.bar(x=centers_mc[0], height=2*stat_err_mc, bottom=np.sum(contents_mc, axis=0) - err_mc, width=width_mc[0], color='grey', hatch='xx', ec='darkgrey', alpha=0.4, lw=0.25, label='MC Uncertainty')
        

    if plot_cfg['include_onbeam'] == True:
        ax.errorbar(centers_onbeam, bins_onbeam, xerr=width_onbeam/2, yerr=stat_err_onbeam, linestyle='none', ecolor='black', label=label_onbeam)

    h, l = ax.get_legend_handles_labels()
    if plot_cfg['include_mc'] == True and plot_cfg['include_mc_err'] == True and plot_cfg['include_onbeam'] == False:
        h_mc, l_mc = h[:-1], l[:-1]
        h_mc = h_mc[::-1]
        h_syst, l_syst = [h[-1]], [l[-1]]
        h, l = h_mc + h_syst, l_mc + l_syst
    elif plot_cfg['include_mc'] == True and plot_cfg['include_mc_err'] == True and plot_cfg['include_onbeam'] == True:
        h_mc, l_mc = h[:-2], l[:-2]
        h_mc = h_mc[::-1]
        h_syst, l_syst = [h[-2]], [l[-2]]
        h_onbeam, l_onbeam = [h[-1]], [l[-1]]
        h, l = h_mc + h_syst + h_onbeam, l_mc + l_syst + l_onbeam

    # Fitting
    fits, fit_labels = [], []
    if plot_cfg['fit_mc'] == True:
        xspace=np.linspace(cfg['variables'][var]['bins'][1], cfg['variables'][var]['bins'][2], 1000)
        popt_mc, pcov_mc = curve_fit(crystalball, centers_mc[0][:], bins_mc[:], p0=(100, 1, 3, 135, 10, -0.01, 0.01))
        y_fit_mc = crystalball(xspace,*popt_mc)
        mc_fit_label = r'$\bf{MC\ Fit}$'+f'\nμ={round(popt_mc[3],1)} ± {round(np.sqrt(np.diag(pcov_mc))[3],1)} [MeV/c$^{2}$]\nσ = {round(popt_mc[4],1)} ± {round(np.sqrt(np.diag(pcov_mc))[4],1)} [MeV/c$^{2}$]'
        mc_fit,=ax.plot(xspace, y_fit_mc, 'b', linewidth=1.5)
        fits.append(mc_fit)
        fit_labels.append(mc_fit_label)
    if plot_cfg['fit_onbeam'] == True:
        xspace=np.linspace(cfg['variables'][var]['bins'][1], cfg['variables'][var]['bins'][2], 1000)
        popt_onbeam, pcov_onbeam = curve_fit(crystalball, centers_onbeam[:], bins_onbeam[:], p0=(135, 1, 3, 135, 10, -0.01, 0.01), maxfev=100000)
        y_fit_onbeam = crystalball(xspace,*popt_onbeam)
        onbeam_fit_label = r'$\bf{Data\ Fit}$'+f'\nμ={round(popt_onbeam[3],1)} ± {round(np.sqrt(np.diag(pcov_onbeam))[3],1)} [MeV/c$^{2}$]\nσ = {round(popt_onbeam[4],1)} \
± {round(np.sqrt(np.diag(pcov_onbeam))[4],1)} [MeV/c$^{2}$]'
        onbeam_fit,=ax.plot(xspace, y_fit_onbeam, 'r', linewidth=1.5)
        fits.append(onbeam_fit)
        fit_labels.append(onbeam_fit_label)
    if plot_cfg['fit_mc'] == True or plot_cfg['fit_onbeam'] == True:
        l1 = ax.legend(fits, fit_labels, fontsize=14, loc='upper right', bbox_to_anchor=(1, 0.68), frameon=False)
        plt.gca().add_artist(l1)

    ax.legend(h,l, ncols=2)
    add_plot_labels(ax, pot_onbeam, adj_y=0.025, title=cfg['variables'][var]['title'])
    ax.set_xlabel(cfg['variables'][var]['xlabel'])
    ax.set_ylabel(cfg['variables'][var]['ylabel'])
    ax.set_ylim(bottom=0)
    plt.show()
        
######################################################################
### Main Analysis
######################################################################
def main(args):

    # Load pi0 config
    pi0_config = toml.load(args.config)

    # Load input file
    rf = uproot.open(args.in_file)
    
    # Load trees
    sel_nu_mc_tree = rf['events/mc/SelectedNu_PhaseCuts']
    sel_nu_mc_df = sel_nu_mc_tree.arrays(library='pd')
    sel_cos_mc_tree = rf['events/mc/SelectedCos_PhaseCuts']
    sel_cos_mc_df = sel_cos_mc_tree.arrays(library='pd')
    sel_mc_df = pd.concat([sel_nu_mc_df, sel_cos_mc_df])
    pot_mc = rf['events/mc/POT'].to_numpy()[0][0]
    sel_onbeam_tree = rf['events/onbeam/SelectedNu_PhaseCuts']
    sel_onbeam_df = sel_onbeam_tree.arrays(library='pd')
    pot_onbeam = rf['events/onbeam/POT'].to_numpy()[0][0]

    # Plotting
    # Config
    plot_config = {'include_mc' : True,
                   'include_mc_err': True,
                   'fit_mc': False,
                   'include_onbeam' : True,
                   'fit_onbeam': False,
                   'normalization' : 'data'}
    #plot_histogram(pi0_config, plot_config, 'pi0_mass', sel_mc_df, pot_mc)
    plot_histogram(pi0_config, plot_config, 'pi0_mass', sel_mc_df, pot_mc, sel_onbeam_df, pot_onbeam)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--in_file', required=True)
    args = parser.parse_args()
    main(args)
