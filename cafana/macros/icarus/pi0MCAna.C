/**
 * @file pi0MCAna.C
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
#include "include/pi0ana/variables_pi0ana_phase.h"
#include "include/pi0ana/cuts_pi0ana_phase.h"
#include "include/pi0ana/variables_pi0ana_nophase.h"
#include "include/pi0ana/cuts_pi0ana_nophase.h"
#include "include/pi0ana/variables_pi0ana_trad.h"
#include "include/pi0ana/cuts_pi0ana_trad.h"
#include "include/spinevar.h"
#include "include/srvar.h"
#include "include/analysis.h"
#include "include/beaminfo.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void pi0MCAna()
{
    //ana::Analysis analysis("pi0ana_bnb_nu_cosmic_mc_cv_v09_89_update");
    ana::Analysis analysis("pi0ana_test");
    
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
    ///////////////////////////////////////
    // Analysis with phase cuts starts here
    ///////////////////////////////////////
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_phase;
    #define CUT cuts::pi0ana_phase::SELCUT
    #define TCUT cuts::neutrino
    vars_selected_nu_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_phase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_phase.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category_topology, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_PhaseCuts", vars_selected_nu_phase, false);

    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_phase;
    #undef TCUT
    #define TCUT cuts::cosmic
    vars_selected_cos_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_phase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_phase.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category_topology, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_PhaseCuts", vars_selected_cos_phase, false);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_phase;
    vars_purity_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(FLASHCUT), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::SELCUT), &cuts::no_cut, &TCUT)});
    analysis.AddTree("Purity_PhaseCuts", vars_purity_phase, false);

    #define SIGCUT cuts::pi0ana_phase::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_phase;
    vars_signal_phase.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"CutType", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::cut_type, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_phase.insert({"IsSignal", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::is_signal, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_phase.insert({"IsData", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::is_not_data, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_phase.insert({"baseline", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pdg", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"cc", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::category, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"interaction_mode", SpineVar<TTYPE,TTYPE>(&vars::neutrino_interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_edep", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"has_single_muon", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::single_muon), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"has_no_charged_pions", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::no_charged_pions), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(FLASHCUT), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::SELCUT), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_PhaseCuts", vars_signal_phase, true);

    // Confusion    
    std::map<std::string, ana::SpillMultiVar> vars_eff_pid_confusion_phase;
    auto true_primary_particle = [](const TTYPEP & p) { return p.is_primary && p.pid >=0; };
    auto true_pid = [](const TTYPEP & p) -> double { return p.pid; };
    auto reco_pid = [](const RTYPEP & p) -> double { return p.pid; };
    vars_eff_pid_confusion_phase.insert({"true_pid", SpineVar<TTYPEP,TTYPEP,TTYPE>(true_pid, true_primary_particle, &SIGCUT)});
    vars_eff_pid_confusion_phase.insert({"reco_pid", SpineVar<RTYPEP,TTYPEP,TTYPE>(reco_pid, true_primary_particle, &SIGCUT)});
    vars_eff_pid_confusion_phase.insert({"reco_primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), true_primary_particle, &SIGCUT)});
    analysis.AddTree("EffPIDConfusion_PhaseCuts", vars_eff_pid_confusion_phase, true);
    
    std::map<std::string, ana::SpillMultiVar> vars_eff_primary_confusion_phase;
    auto true_particle = [](const TTYPEP & p) { return p.pid >=0; };
    vars_eff_primary_confusion_phase.insert({"true_pid", SpineVar<TTYPEP,TTYPEP,TTYPE>(true_pid, true_particle, &SIGCUT)});
    vars_eff_primary_confusion_phase.insert({"reco_pid", SpineVar<RTYPEP,TTYPEP,TTYPE>(reco_pid, true_particle, &SIGCUT)});
    vars_eff_primary_confusion_phase.insert({"true_primary", SpineVar<TTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), true_particle, &SIGCUT)});
    vars_eff_primary_confusion_phase.insert({"reco_primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), true_particle, &SIGCUT)});
    analysis.AddTree("EffPrimaryConfusion_PhaseCuts", vars_eff_primary_confusion_phase, true);

    std::map<std::string, ana::SpillMultiVar> vars_pur_pid_confusion_phase;
    auto reco_primary_particle = [](const RTYPEP & p) { return p.is_primary && p.pid >=0; };
    vars_pur_pid_confusion_phase.insert({"reco_pid", SpineVar<RTYPEP,RTYPEP,RTYPE>(reco_pid, reco_primary_particle, &CUT)});
    vars_pur_pid_confusion_phase.insert({"true_pid", SpineVar<TTYPEP,RTYPEP,RTYPE>(true_pid, reco_primary_particle, &CUT)});
    vars_pur_pid_confusion_phase.insert({"true_primary", SpineVar<TTYPEP,RTYPEP,RTYPE>(WRAP_BOOL(pcuts::is_primary), reco_primary_particle, &CUT)});
    analysis.AddTree("PurPIDConfusion_PhaseCuts", vars_pur_pid_confusion_phase, true);
    
    //////////////////////////////////////////
    // Analysis without phase cuts starts here
    //////////////////////////////////////////
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_nophase;
    #undef CUT
    #define CUT cuts::pi0ana_nophase::SELCUT
    #undef TCUT
    #define TCUT cuts::neutrino
    vars_selected_nu_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nophase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nophase.insert({"IsData", SrVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)}); // GUNDAM 
    vars_selected_nu_nophase.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category_topology, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_NoPhaseCuts", vars_selected_nu_nophase, true);

    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_nophase;
    #undef TCUT
    #define TCUT cuts::cosmic
    vars_selected_cos_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nophase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nophase.insert({"IsData", SrVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nophase.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category_topology, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_NoPhaseCuts", vars_selected_cos_nophase, true);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_nophase;
    vars_purity_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(FLASHCUT), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::SELCUT), &cuts::no_cut, &TCUT)});
    analysis.AddTree("Purity_NoPhaseCuts", vars_purity_nophase, false);

    #undef SIGCUT
    #define SIGCUT cuts::pi0ana_nophase::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_nophase;
    vars_signal_nophase.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"baseline", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pdg", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"cc", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::category, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"interaction_mode", SpineVar<TTYPE,TTYPE>(&vars::neutrino_interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_edep", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"has_single_muon", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::single_muon), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"has_no_charged_pions", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::no_charged_pions), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(FLASHCUT), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::SELCUT), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_NoPhaseCuts", vars_signal_nophase, true);    

    /////////////////////////////////////////////
    // Analysis with traditional cuts starts here
    /////////////////////////////////////////////
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_trad;
    #undef CUT
    #define CUT cuts::pi0ana_trad::SELCUT
    #undef TCUT
    #define TCUT cuts::neutrino
    vars_selected_nu_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_trad.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_trad.insert({"IsData", SrVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_momentum_mag, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_trad.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category_topology, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_TradCuts", vars_selected_nu_trad, true);

    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_trad;
    #undef TCUT
    #define TCUT cuts::cosmic
    vars_selected_cos_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"IsData", SrVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_momentum_mag, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"baseline", SpineVar<TTYPE,RTYPE>(&vars::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category_topology, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"interaction_mode", SpineVar<TTYPE,RTYPE>(&vars::neutrino_interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_TradCuts", vars_selected_cos_trad, true);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_trad;
    vars_purity_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(FLASHCUT), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::SELCUT), &cuts::no_cut, &TCUT)});
    analysis.AddTree("Purity_TradCuts", vars_purity_trad, false);

    #undef SIGCUT
    #define SIGCUT cuts::pi0ana_trad::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_trad;
    vars_signal_trad.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"baseline", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pdg", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"cc", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::category, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"interaction_mode", SpineVar<TTYPE,TTYPE>(&vars::neutrino_interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_edep", SpineVar<TTYPE,TTYPE>(&vars::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"has_single_muon", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::single_muon), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"has_no_charged_pions", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::no_charged_pions), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(FLASHCUT), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::SELCUT), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_TradCuts", vars_signal_trad, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
