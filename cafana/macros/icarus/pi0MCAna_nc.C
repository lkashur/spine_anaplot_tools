/**
 * @file pi0MCAna_nc.C
 * @brief The main analysis macro for the pi0ana analysis on ICARUS Monte
 * Carlo simulation.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author lkashur@colostate.edu
*/
#include "include/variables.h"
#include "include/cuts.h"
#include "include/pi0ana/variables_pi0ana_nc.h"
#include "include/pi0ana/cuts_pi0ana_nc.h"
#include "include/spinevar.h"
#include "include/srvar.h"
#include "include/analysis.h"
#include "include/beaminfo_nc.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void pi0MCAna_nc()
{
    //ana::Analysis analysis("pi0ana_bnb_nu_cosmic_mc_cv_v09_89_update");
    ana::Analysis analysis("pi0ana_nc_test");
    
    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/fall2024/collonly_v2b/flat/*.root"); // BNB Mini Sample
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/spineprod/mcsim/nominal/flat/*.root"); // Nominal BNB Nu + Cosmic (v09_89)
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/lkashur/v09_89_01_01p02_numi_nu_cosmic_mc/flat/*.root"); // NuMI

    //#define IS_NUMI 1 // UNCOMMENT IF PROCESSING NUMI MC
    
    analysis.AddLoader("mc", &mc, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    /////////////////////////////////////////////
    // NC analysis starts here
    /////////////////////////////////////////////
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_nc;
    #undef CUT
    #define CUT cuts::pi0ana_nc::SELCUT
    #undef TCUT
    #define TCUT cuts::neutrino
    vars_selected_nu_nc.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nc.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nc::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nc.insert({"IsData", SrVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_momentum_mag, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nc.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nc::category, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nc::category_topology, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_nu_nc.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_NCCuts", vars_selected_nu_nc, true);

    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_nc;
    #undef TCUT
    #define TCUT cuts::cosmic
    vars_selected_cos_nc.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nc.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nc::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nc.insert({"IsData", SrVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_momentum_mag, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nc.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nc::category, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nc::category_topology, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nc::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_cos_nc.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_NCCuts", vars_selected_cos_nc, true);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_nc;
    vars_purity_nc.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_nc.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nc::category, &cuts::no_cut, &TCUT)});
    vars_purity_nc.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nc::topological_0mu2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nc::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(FLASHCUT), &cuts::no_cut, &TCUT)});
    vars_purity_nc.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nc::SELCUT), &cuts::no_cut, &TCUT)});
    analysis.AddTree("Purity_NCCuts", vars_purity_nc, false);

    #undef SIGCUT
    #define SIGCUT cuts::pi0ana_nc::signal_0mu1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_nc;
    vars_signal_nc.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"baseline", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pdg", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"cc", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::category, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"interaction_mode", SpineVar<TTYPE,TTYPE>(&vars::neutrino_interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"true_edep", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nc::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nc::topological_0mu2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nc::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"has_single_muon", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nc::single_muon), &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"has_no_charged_pions", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nc::no_charged_pions), &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(FLASHCUT), &SIGCUT, &SIGCUT)});
    vars_signal_nc.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nc::SELCUT), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_NCCuts", vars_signal_nc, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
