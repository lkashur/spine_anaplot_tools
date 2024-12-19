/**
 * @file pi02024.C
 * @brief The main analysis macro for the pi02024 analysis on ICARUS Monte
 * Carlo simulation.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author lkashur@colostate.edu
*/
#include "include/variables.h"
#include "include/pi02024/variables_pi02024.h"
#include "include/cuts.h"
#include "include/pi02024/cuts_pi02024.h"
#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void pi02024()
{
    ana::Analysis analysis("pi02024_rev2_icarus_v4");

    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/fall2024/collonly_v2b/flat/*.root");
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/fall2024/nominal/flat/*.root");
    analysis.AddLoader("mc", &mc, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::pi02024::all_1mu0pi2gamma_cut
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu;
    vars_selected_nu.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi02024::category, &CUT, &TCUT)});
    vars_selected_nu.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi02024::category_topology, &CUT, &TCUT)});
    vars_selected_nu.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_nu.insert({"muon_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::pi02024::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi02024::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_nu.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_nu.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});

    analysis.AddTree("selectedNu", vars_selected_nu, false);

    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define SIGCUT cuts::pi02024::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal;
    vars_signal.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"baseline", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pdg", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"cc", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi02024::category, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi02024::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"interaction_mode", SpineVar<TTYPE,TTYPE>(&vars::neutrino_interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi02024::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi02024::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"flash_cut_bnb", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut_bnb), &SIGCUT, &SIGCUT)});

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
