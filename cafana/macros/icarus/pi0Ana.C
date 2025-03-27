/**
 * @file pi0MCAna.C
 * @brief The main analysis macro for the ICARUS numu CC pi0 selection..
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author lkashur@colostate.edu
*/
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PIDFUNC pvars::pid
//#define PIDFUNC pvars::custom_pid
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false
#define WRITE_PURITY_TREES false

#include "include/mctruth.h"
#include "include/variables.h"
#include "include/cuts.h"
#include "include/pi0ana/variables_pi0ana_phase.h"
#include "include/pi0ana/cuts_pi0ana_phase.h"
#include "include/pi0ana/variables_pi0ana_nophase.h"
#include "include/pi0ana/cuts_pi0ana_nophase.h"
#include "include/pi0ana/variables_pi0ana_trad.h"
#include "include/pi0ana/cuts_pi0ana_trad.h"

#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void pi0Ana()
{

    ana::Analysis analysis("offbeam_27_mar_2025");

    //ana::SpectrumLoader onbeam("/pnfs/icarus/persistent/users/mueller/production/data/onbeam/flat/input*.flat.root");
    //analysis.AddLoader("onbeam", &onbeam, false);

    ana::SpectrumLoader offbeam("/pnfs/icarus/persistent/users/mueller/production/data/offbeam/flat/input*.flat.root"); 
    analysis.AddLoader("offbeam", &offbeam, false);

    //ana::SpectrumLoader run9435("/pnfs/icarus/persistent/users/mueller/production/data/onbeam_run9435/flat/input.flat.root");
    //analysis.AddLoader("run9435", &run9435, false);
                         
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/production/simulation/nominal/flat/input*.flat.root");
    //analysis.AddLoader("mc", &mc, true);




    ///////////////////////////////////////////
    /// Monte Carlo (CV)
    ///////////////////////////////////////////
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/lkashur/bnb_nu_cosmic_cv_merged/flatcaf*.root"); // BNB
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/production/simulation/nominal/flat/input*.flat.root"); // updated CV"
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/lkashur/NuMI_nu_v09_89_01_01p03/flat/*.root"); // NuMI
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/lkashur/NuMI_nu_v09_89_01_01p03/flat/prodcorsika_genie_protononly_icarus_numi_gen_g4_71567695_0-163b2f1e-c022-4045-a546-4db951aa37fd.spine.caf.flat.root");
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/production/simulation/nominal/flat/input*.root"); // BNB nu + cosmic CV, upstream calibration fixed

    //ana::SpectrumLoader var01("/pnfs/icarus/persistent/users/mueller/production/simulation/var01/flat/input*.root");
    
    //analysis.AddLoader("mc", &mc, true);
    //analysis.AddLoader("var01", &var01, true);

    ///////////////////////////////////////////
    /// Monte Carlo (var0)
    ///////////////////////////////////////////

    ///////////////////////////////////////////
    /// Monte Carlo (var1)
    ///////////////////////////////////////////

    /**
     * @brief Add variabls for selected interactions (in-phase) to the analysis.
     * @details This adds a set of variables to the analysis by creating a map
     * of variable names and SpillMultiVars that provide the functionality to
     * create the variables.  These names are used in the TTree that is created
     * by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::pi0ana_phase::all_1mu0pi2gamma_cut
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_phase;
    vars_selected_nu_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_phase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_phase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category_topology, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"category_topology_simple", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category_topology_simple, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"visible_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::visible_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_PhaseCuts", vars_selected_nu_phase, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_phase;
    vars_selected_cos_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_phase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_phase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category_topology, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"category_topology_simple", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category_topology_simple, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"visible_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_phase::visible_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_PhaseCuts", vars_selected_cos_phase, false);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_phase;
    vars_purity_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_phase::category, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"one_muon_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::one_muon_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"zero_charged_pions_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::zero_charged_pions_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"two_or_three_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::two_or_three_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_phase::all_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    if constexpr(WRITE_PURITY_TREES)
		  analysis.AddTree("Purity_PhaseCuts", vars_purity_phase, false);
    //analysis.AddTree("Purity_PhaseCuts", vars_purity_phase, true);

    #define SIGCUT cuts::pi0ana_phase::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_phase;
    vars_signal_phase.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"CutType", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::cut_type, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_phase.insert({"IsSignal", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::is_signal, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_phase.insert({"IsData", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::is_not_data, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_phase.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::category, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_photons_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_phase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"one_muon_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::one_muon_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"zero_charged_pions_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::zero_charged_pions_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"two_or_three_photons_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::two_or_three_photons_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_vertex_x", SpineVar<TTYPE,TTYPE>(&vars::vertex_x, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_vertex_y", SpineVar<TTYPE,TTYPE>(&vars::vertex_y, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_vertex_z", SpineVar<TTYPE,TTYPE>(&vars::vertex_z, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_phase::all_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_PhaseCuts", vars_signal_phase, true);

    /**
     * @brief Add variabls for selected interactions (no-phase) to the analysis.
     * @details This adds a set of variables to the analysis by creating a map
     * of variable names and SpillMultiVars that provide the functionality to
     * create the variables.  These names are used in the TTree that is created
     * by the Tree class to store the results of the analysis.
     */
    #undef CUT
    #define CUT cuts::pi0ana_nophase::all_1mu0pi2gamma_cut
    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_nophase;
    vars_selected_nu_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nophase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nophase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nophase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category_topology, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_NoPhaseCuts", vars_selected_nu_nophase, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_nophase;
    vars_selected_cos_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nophase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nophase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nophase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category_topology, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_nophase::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_NoPhaseCuts", vars_selected_cos_nophase, false);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_nophase;
    vars_purity_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_nophase::category, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"one_muon_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::one_muon_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"zero_charged_pions_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::zero_charged_pions_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"two_or_three_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::two_or_three_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::all_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    if constexpr(WRITE_PURITY_TREES)
		  analysis.AddTree("Purity_NoPhaseCuts", vars_purity_nophase, false);

    #undef SIGCUT
    #define SIGCUT cuts::pi0ana_nophase::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_nophase;
    vars_signal_nophase.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"CutType", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::cut_type, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_nophase.insert({"IsSignal", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::is_signal, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_nophase.insert({"IsData", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::is_not_data, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_nophase.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::category, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_photons_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_nophase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"one_muon_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::one_muon_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"zero_charged_pions_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::zero_charged_pions_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"two_or_three_photons_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::two_or_three_photons_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_vertex_x", SpineVar<TTYPE,TTYPE>(&vars::vertex_x, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_vertex_y", SpineVar<TTYPE,TTYPE>(&vars::vertex_y, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_vertex_z", SpineVar<TTYPE,TTYPE>(&vars::vertex_z, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_nophase::all_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_NoPhaseCuts", vars_signal_nophase, true);

    /**
     * @brief Add variabls for selected interactions (traditional) to the analysis.
     * @details This adds a set of variables to the analysis by creating a map
     * of variable names and SpillMultiVars that provide the functionality to
     * create the variables.  These names are used in the TTree that is created
     * by the Tree class to store the results of the analysis.
     */
    #undef CUT
    #define CUT cuts::pi0ana_trad::all_1mu0pi2gamma_cut
    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_trad;
    vars_selected_nu_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_trad.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_trad.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::is_not_data, &CUT, &TCUT)}); // GUNDAM 
    vars_selected_nu_trad.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category_topology, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"category_topology_simple", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category_topology_simple, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_TradCuts", vars_selected_nu_trad, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_trad;
    vars_selected_cos_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category_topology, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"category_topology_simple", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category_topology_simple, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::pi0ana_trad::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_TradCuts", vars_selected_cos_trad, false);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_trad;
    vars_purity_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::pi0ana_trad::category, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"one_muon_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::one_muon_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"zero_charged_pions_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::zero_charged_pions_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"two_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::two_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"two_or_three_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::two_or_three_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::pi0ana_trad::all_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    if constexpr(WRITE_PURITY_TREES)
		  analysis.AddTree("Purity_TradCuts", vars_purity_trad, false);

    #undef SIGCUT
    #define SIGCUT cuts::pi0ana_trad::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_trad;
    vars_signal_trad.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"CutType", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::cut_type, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_trad.insert({"IsSignal", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::is_signal, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_trad.insert({"IsData", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::is_not_data, &SIGCUT, &SIGCUT)}); // GUNDAM 
    vars_signal_trad.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::category, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_photons_costheta", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::pi0ana_trad::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"one_muon_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::one_muon_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"zero_charged_pions_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::zero_charged_pions_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"two_or_three_photons_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::two_or_three_photons_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_vertex_x", SpineVar<TTYPE,TTYPE>(&vars::vertex_x, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_vertex_y", SpineVar<TTYPE,TTYPE>(&vars::vertex_y, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_vertex_z", SpineVar<TTYPE,TTYPE>(&vars::vertex_z, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::pi0ana_trad::all_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
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
