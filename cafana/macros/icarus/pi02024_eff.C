/**
 * @file pi02024_eff.C
 * @brief The main analysis macro for the pi02024 analysis on ICARUS Monte
 * Carlo simulation.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author lkashur@colostate.edu
*/
#include "include/variables.h"
#include "include/pi02024/variables_pi02024_eff.h"
#include "include/cuts.h"
#include "include/pi02024/cuts_pi02024_eff.h"
#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void pi02024_eff()
{
    ana::Analysis analysis("pi02024_eff");

    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/fall2024/collonly_v2b/flat/*.root");
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/fall2024/nominal/flat/*.root");
    analysis.AddLoader("mc", &mc, true);

    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_purity_nu;
    vars_purity_nu.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi02024_eff::category, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi02024_eff::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"flash_cut_bnb", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut_bnb), &cuts::no_cut, &TCUT)});
    
    analysis.AddTree("purityNu", vars_purity_nu, false);

    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define SIGCUT cuts::pi02024_eff::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal;
    vars_signal.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"baseline", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pdg", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"cc", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi02024_eff::category, &SIGCUT, &SIGCUT)});
    //vars_signal.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi02024::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"interaction_mode", SpineVar<TTYPE,TTYPE>(&vars::neutrino_interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_edep", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi02024_eff::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi02024_eff::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi02024_eff::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi02024_eff::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi02024_eff::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi02024_eff::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"topological_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi02024_eff::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"flash_cut_bnb", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut_bnb), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi02024_eff::all_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});

    analysis.AddTree("signal", vars_signal, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
