/**
 * @file pi0AnaNC.C
 * @brief The main analysis macro for the ICARUS numu NC pi0 selection..
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

//#include "include/pi0ana/variables_pi0ana_nc_inc.h"
//#include "include/pi0ana/cuts_pi0ana_nc_inc.h"
//#include "include/pi0ana/variables_pi0ana_nc_inc.h"
//#define SEL pi0ana_nc_inc

#include "include/pi0ana/variables_pi0ana_nc_0p.h"
#include "include/pi0ana/cuts_pi0ana_nc_0p.h"
#include "include/pi0ana/variables_pi0ana_nc_0p.h"
#define SEL pi0ana_nc_0p

//#include "include/pi0ana/variables_pi0ana_nc_g1p.h"
//#include "include/pi0ana/cuts_pi0ana_nc_g1p.h"
//#include "include/pi0ana/variables_pi0ana_nc_g1p.h"
//#define SEL pi0ana_nc_g1p

#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void pi0AnaNC()
{

    ana::Analysis analysis("test_nc_0p_28_mar_2025");

    //ana::SpectrumLoader onbeam("/pnfs/icarus/persistent/users/mueller/production/data/onbeam/flat/input*.flat.root");
    //analysis.AddLoader("onbeam", &onbeam, false);

    //ana::SpectrumLoader offbeam("/pnfs/icarus/persistent/users/mueller/production/data/offbeam/flat/input*.flat.root"); 
    //analysis.AddLoader("offbeam", &offbeam, false);

    //ana::SpectrumLoader run9435("/pnfs/icarus/persistent/users/mueller/production/data/onbeam_run9435/flat/input.flat.root");
    //analysis.AddLoader("run9435", &run9435, false);
                         
    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/production/simulation/nominal/flat/input*.flat.root");
    analysis.AddLoader("mc", &mc, true);

    /**
     * @brief Add variabls for selected interactions (in-phase) to the analysis.
     * @details This adds a set of variables to the analysis by creating a map
     * of variable names and SpillMultiVars that provide the functionality to
     * create the variables.  These names are used in the TTree that is created
     * by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::SEL::all_0mu0pi2gamma_cut
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_nc_inc;
    vars_selected_nu_nc_inc.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::SEL::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nc_inc.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::SEL::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nc_inc.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::SEL::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_nc_inc.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::SEL::category, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::SEL::category_topology, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_mass, &CUT, &TCUT)});
    //vars_selected_nu_nc_inc.insert({"visible_energy", SpineVar<RTYPE,RTYPE>(&vars::SEL::visible_energy, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_nc_inc.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_NC_Cuts", vars_selected_nu_nc_inc, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_nc_inc;
    vars_selected_cos_nc_inc.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::SEL::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nc_inc.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::SEL::is_signal, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nc_inc.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::SEL::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_nc_inc.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::SEL::category, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::SEL::category_topology, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::SEL::pi0_mass, &CUT, &TCUT)});
    //vars_selected_cos_nc_inc.insert({"visible_energy", SpineVar<RTYPE,RTYPE>(&vars::SEL::visible_energy, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_nc_inc.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_NC_Cuts", vars_selected_cos_nc_inc, false);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_nc_inc;
    vars_purity_nc_inc.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::SEL::category, &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"one_muon_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::SEL::one_muon_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"zero_charged_pions_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::SEL::zero_charged_pions_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"two_or_three_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::SEL::two_or_three_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::SEL::topological_0mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::SEL::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::SEL::all_0mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_nc_inc.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    if constexpr(WRITE_PURITY_TREES)
		  analysis.AddTree("Purity_NC_Cuts", vars_purity_nc_inc, false);

    #define SIGCUT cuts::SEL::signal_0mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_nc_inc;
    vars_signal_nc_inc.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"CutType", SpineVar<TTYPE,TTYPE>(&vars::SEL::cut_type, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_nc_inc.insert({"IsSignal", SpineVar<TTYPE,TTYPE>(&vars::SEL::is_signal, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_nc_inc.insert({"IsData", SpineVar<TTYPE,TTYPE>(&vars::SEL::is_not_data, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_nc_inc.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::SEL::category, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::SEL::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::SEL::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::SEL::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::SEL::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::SEL::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::SEL::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::SEL::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_photons_costheta", SpineVar<TTYPE,TTYPE>(&vars::SEL::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::SEL::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"one_muon_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::SEL::one_muon_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"zero_charged_pions_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::SEL::zero_charged_pions_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"two_or_three_photons_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::SEL::two_or_three_photons_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::SEL::topological_0mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::SEL::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"true_vertex_x", SpineVar<TTYPE,TTYPE>(&vars::vertex_x, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"true_vertex_y", SpineVar<TTYPE,TTYPE>(&vars::vertex_y, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"true_vertex_z", SpineVar<TTYPE,TTYPE>(&vars::vertex_z, &SIGCUT, &SIGCUT)});
    vars_signal_nc_inc.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::SEL::all_0mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_NC_Cuts", vars_signal_nc_inc, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
