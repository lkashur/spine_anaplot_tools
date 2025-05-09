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
    rf = uproot.open(args.in_file)
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
    ########################################

    ##############################################
    ### Register variables from MC sample
    mc_sample = samples['mc']
    for v in variables.values():
        mc_sample.register_variable(v, categories)
    ##############################################

    ##############################################
    ### Get covariance matrix from each knob
    for k,v in mc_sample._systematics.items():
        v.process(mc_sample)
    ##############################################

    ##############################################
    ### Systematics
    ##############################################
    if args.syst == 'xsec':
        regxp = re.compile('\\bGENIEReWeight_SBN_v1_multisim_\\w+\\b')
    elif args.syst == 'flux':
        regxp = re.compile('\\b\\w*_Flux\\b')

    # All knobs, individually
    systs = [syst for syst in mc_sample._systematics.values() if regxp.match(syst._name)]
    

    #############################
    ### Cross Section Systematics
    #############################
    # All knobs, individually
    regxp = re.compile('\\bGENIEReWeight_SBN_v1_multisim_\\w+\\b')
    xsec_systs = [syst for syst in mc_sample._systematics.values() if regxp.match(syst._name)]
    #for syst in xsec_systs:
    #    for k,v in syst._covariances.items():
    #        print(k)
    #        print(v)

    # Combine knobs
    xsec_syst_total = Systematic.combine(xsec_systs, 'xsec', 'Cross Section Systematics')
    mc_sample._systematics[xsec_syst_total._name] = xsec_syst_total

    # Get x-sec covariance matrix for a single variable
    xsec_var_cov = mc_sample._systematics['xsec']._covariances[f'xsec_{args.var}']
 
    # Convert to TMatrixTSym
    xsec_var_cov_tmatrix = ndarray_to_tmatrixtsym(xsec_var_cov)

    ####################
    ### Flux Systematics
    ####################
    
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
    parser.add_argument('--in_file', required=False)
    parser.add_argument('--syst', required=False)
    parser.add_argument('--var', required=False)
    args = parser.parse_args()
    main(args)
