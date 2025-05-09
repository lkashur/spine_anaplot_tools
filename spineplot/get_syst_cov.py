# imports
import os
import toml
import uproot
import argparse
import re
from matplotlib import pyplot as plt
import ROOT

from sample import Sample
from figure import SpineFigure, SimpleFigure
from spectra1d import SpineSpectra1D
from spectra2d import SpineSpectra2D
from efficiency import SpineEfficiency
from roc import ROCCurve
from ternary import Ternary
from style import Style
from variable import Variable
from systematic import Systematic

def main(args):

    ##############################
    ### Load analysis file
    rf = uproot.open(args.input)
    ##############################

    ###############################
    ### Load config
    config = toml.load(args.config)
    ###############################

    #######################################
    ### Load samples and systematic recipes
    samples = {name: Sample(name, rf, config['analysis']['category_branch'], **config['samples'][name]) for name in config['samples']}
    variables = {name: Variable(name, **config['variables'][name]) for name in config['variables']}
    categories = dict()
    for ci, cat in enumerate(config['analysis']['category_assignment']):
        categories.update({c : config['analysis']['category_labels'][ci] for c in cat})
    recipes = config['systematic_recipe']


    # Test
    for s in samples.values():
        for v in variables.values():
            s.register_variable(v, categories)
        s.process_systematics(recipes)
        s.set_weight(target=config['analysis']['ordinate_sample'])
        systs = [syst for syst in s._systematics.values()]
        print(systs)
        
    print('hello')
    ########################################

    ##############################################
    ### Register variables from MC sample
    mc_nu_sample = samples['mc_nu']
    mc_cos_sample = samples['mc_cos']
    mc_offbeam_sample = samples['offbeam']
    for v in variables.values():
        mc_nu_sample.register_variable(v, categories)
        mc_cos_sample.register_variable(v, categories)
        mc_offbeam_sample.register_variable(v, categories)
    ##############################################

    ##############################################
    ### Get covariance matrix from each parameter
    for k,v in mc_nu_sample._systematics.items():
        v.process(mc_nu_sample, mc_nu_sample._presel_mask)

    for k,v in mc_cos_sample._systematics.items():
        v.process(mc_cos_sample, None)

    for k,v in mc_offbeam_sample._systematics.items():
        v.process(mc_offbeam_sample, None)
    ##############################################

    ##############################################
    ### Systematics
    ##############################################
    if args.syst == 'xsec':
        #regxp = re.compile('\\bGENIEReWeight_SBN_v1_multisim_\\w+\\b')
        regxp = re.compile('GENIEReWeight_SBN_v1_multisim')
    elif args.syst == 'flux':
        #regxp = re.compile('\\b\\w*_Flux\\b')
        regxp = re.compile('_Flux')
    elif args.syst == 'det':
        regxp = re.compile('var(0|1|2|3|4|5|6|7|8|9|10)+')
    elif args.syst == 'stat':
        regxp = re.compile('statistical')

    # All knobs, individually
    systs = [syst for syst in mc_nu_sample._systematics.values() if regxp.search(syst._name)]

    for syst in systs:
        print(syst)
    
    # Generic
    chosen_systs = [syst for syst in mc_nu_sample._systematics.values() if regxp.search(syst._name)]
    total_chosen_syst = Systematic.combine(chosen_systs, 'total_chosen_syst', None)
    chosen_fractional_uncertainties = []
    for syst in chosen_systs:
        chosen_fractional_uncertainties.append({'name':syst._name, 'std':round(100*syst._std, 1)})
    #print(chosen_fractional_uncertainties)
    #print(total_chosen_syst)
    
    #############################
    ### Cross Section Systematics
    #############################
    '''
    # All knobs, individually
    xsec_systs = [syst for syst in mc_sample._systematics.values() if regxp.search(syst._name)]
    total_xsec_syst = Systematic.combine(xsec_systs, 'total_xsec_syst', None)
    print(total_xsec_syst)
    xsec_fractional_uncertainties = []
    for syst in xsec_systs:
        
        #for k,v in syst._covariances.items():
        #    print(k)
        #    print(v)
        
        xsec_fractional_uncertainties.append({'name':syst._name, 'std':round(100*syst._std, 1)})
    
    #print(xsec_fractional_uncertainties)
    xsec_fractional_uncertainties = sorted(xsec_fractional_uncertainties, key=lambda x: x['std'])
    for i in xsec_fractional_uncertainties:
        print(i)

    # Combine knobs
    xsec_syst_total = Systematic.combine(xsec_systs, 'xsec', 'Cross Section Systematics')
    mc_sample._systematics[xsec_syst_total._name] = xsec_syst_total

    # Get x-sec covariance matrix for a single variable
    xsec_var_cov = mc_sample._systematics['xsec']._covariances[f'xsec_{args.var}']
 
    # Convert to TMatrixTSym
    xsec_var_cov_tmatrix = ndarray_to_tmatrixtsym(xsec_var_cov)
    '''

    ####################
    ### Flux Systematics
    ####################
    '''
    flux_systs = [syst for syst in mc_sample._systematics.values() if regxp.search(syst._name)]
    total_flux_syst = Systematic.combine(flux_systs, 'total_flux_syst', None)
    print(total_flux_syst)
    flux_fractional_uncertainties = []
    for syst in flux_systs:
        flux_fractional_uncertainties.append({'name':syst._name, 'std':round(100*syst._std, 1)})
    flux_fractional_uncertainties = sorted(flux_fractional_uncertainties, key=lambda x: x['std'])
    for i in flux_fractional_uncertainties:
        print(i) 
    '''

    ####################
    ### Output TFile
    ####################
    #outf = ROOT.TFile(f'covmat_{args.var}.root', 'RECREATE')
    #outf.WriteObject(xsec_var_cov_tmatrix, f'xsec_{args.var}_cov')
    #outf.Close()
    
    

# Function to convert ndarray to TMatrixTSym<double>
def ndarray_to_tmatrixtsym(np_array):
    size = np_array.shape[0]
    tmatrix = ROOT.TMatrixDSym(size)
    
    for i in range(size):
        for j in range(i, size):
            tmatrix[i][j] = np_array[i, j]
            tmatrix[j][i] = np_array[i, j]
    return tmatrix

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=False)
    parser.add_argument('--input', required=False)
    parser.add_argument('--syst', required=False)
    parser.add_argument('--var', required=False)
    args = parser.parse_args()
    main(args)
