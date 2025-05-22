/**
 * @file ccpi0AnaMC.C
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
#define CALOKEFUNC pvars::calo_ke
//#define CALOKEFUNC pvars::custom_calo_ke
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false
#define WRITE_PURITY_TREES false

#include "include/mctruth.h"
#include "include/variables.h"
#include "include/cuts.h"

#include "include/pi0ana/variables_ccpi0ana_phase.h"
#include "include/pi0ana/cuts_ccpi0ana_phase.h"
//#include "include/pi0ana/variables_ccpi0ana_nophase.h"
//#include "include/pi0ana/cuts_ccpi0ana_nophase.h"
//#include "include/pi0ana/variables_ccpi0ana_trad.h"
//#include "include/pi0ana/cuts_ccpi0ana_trad.h"

#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

/**
 * Usage:
 * ./ccpi0Ana <input file pattern> <output file name> <directory> <data_or_sim>
 * <input file pattern> is the path to input flat CAFs
 * <output file name> is name of saved output ROOT file
 * <directory> refers to output directory to save TTrees (mc, offbeam, onbeam)
 * "data_or_sim" refers to type of input
 */
int main(int argc, char ** argv)
{
    // Output filename
    ana::Analysis analysis(argv[2]);

    // Input
    ana::SpectrumLoader sl(argv[1] + std::string("*flat.root"));

    // Configure loader
    analysis.AddLoader(argv[3], &sl, std::string(argv[4]) == "sim" ? true : false);

    /**
     * @brief Add variabls for selected interactions (in-phase) to the analysis.
     * @details This adds a set of variables to the analysis by creating a map
     * of variable names and SpillMultiVars that provide the functionality to
     * create the variables.  These names are used in the TTree that is created
     * by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::ccpi0ana_phase::all_cut
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_phase;
    vars_selected_nu_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// GUNDAM VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////////
    vars_selected_nu_phase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::cut_type, &CUT, &TCUT)});
    if(std::string(argv[4]) == "sim")
    {
        vars_selected_nu_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::is_not_data, &CUT, &TCUT)});
        vars_selected_nu_phase.insert({"IsNu", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::is_nu, &CUT, &TCUT)});
    }
    else
    {
        vars_selected_nu_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::is_data, &CUT, &TCUT)});
	vars_selected_nu_phase.insert({"IsNu", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::is_not_nu, &CUT, &TCUT)});
    }
    /// END GUNDAM VARIABLES /////////////////////////////////////////////////////////////////////////////////////////////

    vars_selected_nu_phase.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::category, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_muon_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_muon_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_leading_photon_cosphi", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_leading_photon_cosphi", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_leading_photon_ip", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_ip, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_leading_photon_ip", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_ip, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_subleading_photon_cosphi", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_subleading_photon_cosphi", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_subleading_photon_ip", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_subleading_photon_ip", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_photons_avg_ip", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_photons_avg_ip, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_photons_avg_ip", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_photons_avg_ip, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_pi0_mass", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_visible_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::visible_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_visible_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::visible_energy, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_nu_phase.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_PhaseCuts", vars_selected_nu_phase, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_phase;
    vars_selected_cos_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// GUNDAM VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////////
    vars_selected_cos_phase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::cut_type, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"IsNu", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::is_not_nu, &CUT, &TCUT)});
    if(std::string(argv[4]) == "sim")
      vars_selected_cos_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::is_not_data, &CUT, &TCUT)});
    else
      vars_selected_cos_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::is_data, &CUT, &TCUT)});
    /// END GUNDAM VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////
    
    vars_selected_cos_phase.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::category, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_muon_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_muon_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_leading_photon_cosphi", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_leading_photon_cosphi", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_leading_photon_ip", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_ip, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_leading_photon_ip", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_ip, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_subleading_photon_cosphi", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_subleading_photon_cosphi", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_subleading_photon_ip", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_ip, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_subleading_photon_ip", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_ip, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_photons_avg_ip", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_photons_avg_ip, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_photons_avg_ip", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_photons_avg_ip, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_pi0_mass", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_visible_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::visible_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_visible_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::visible_energy, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_cos_phase.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_PhaseCuts", vars_selected_cos_phase, false);

    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_phase;
    vars_purity_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::category, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"one_muon_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::one_muon_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"zero_charged_pions_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::zero_charged_pions_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"two_or_three_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::two_or_three_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"base_topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::base_topology_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"leading_shower_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::leading_shower_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::all_cut), &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_phase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    if constexpr(WRITE_PURITY_TREES)
		  analysis.AddTree("Purity_PhaseCuts", vars_purity_phase, false);

    #define SIGCUT cuts::ccpi0ana_phase::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_phase;
    vars_signal_phase.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// GUNDAM VARIABLES /////////////////////////////////////////////////////////////////////////////////////////////
    vars_signal_phase.insert({"CutType", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::cut_type, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"IsData", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::is_not_data, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"IsNu", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::is_nu, &SIGCUT, &SIGCUT)});
    /// END GUNDAM VARIABLES /////////////////////////////////////////////////////////////////////////////////////////
    
    vars_signal_phase.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::category, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"reco_muon_momentum_mag", SpineVar<RTYPE,TTYPE>(&vars::ccpi0ana_phase::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"reco_muon_beam_costheta", SpineVar<RTYPE,TTYPE>(&vars::ccpi0ana_phase::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_leading_photon_cosphi", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_cosphi, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_leading_photon_ip", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_leading_photon_ip, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_subleading_photon_cosphi", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_cosphi, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_subleading_photon_ip", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_subleading_photon_ip, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"reco_pi0_momentum_mag", SpineVar<RTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"reco_pi0_beam_costheta", SpineVar<RTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_photons_avg_ip", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_phase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"base_topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::base_topology_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"leading_shower_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::leading_shower_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"one_muon_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::one_muon_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"zero_charged_pions_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::zero_charged_pions_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"two_or_three_photons_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::two_or_three_photons_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_vertex_x", SpineVar<TTYPE,TTYPE>(&vars::vertex_x, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_vertex_y", SpineVar<TTYPE,TTYPE>(&vars::vertex_y, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"true_vertex_z", SpineVar<TTYPE,TTYPE>(&vars::vertex_z, &SIGCUT, &SIGCUT)});
    vars_signal_phase.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_phase::all_cut), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_PhaseCuts", vars_signal_phase, true);

    /**
     * @brief Add variabls for selected interactions (no-phase) to the analysis.
     * @details This adds a set of variables to the analysis by creating a map
     * of variable names and SpillMultiVars that provide the functionality to
     * create the variables.  These names are used in the TTree that is created
     * by the Tree class to store the results of the analysis.
     */
    /*
    #undef CUT
    #define CUT cuts::ccpi0ana_nophase::all_1mu0pi2gamma_cut
    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_nophase;
    vars_selected_nu_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    /// GUNDAM VARIABLES /////////////////////////////////////////////////////////////////////////////////////////////////////
    vars_selected_nu_nophase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::cut_type, &CUT, &TCUT)});
    if(std::string(argv[4]) == "sim")
    {
        vars_selected_nu_nophase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::is_not_data, &CUT, &TCUT)});
        vars_selected_nu_nophase.insert({"IsNu", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::is_nu, &CUT, &TCUT)});
    }
    else
    {
        vars_selected_nu_nophase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::is_data, &CUT, &TCUT)});
        vars_selected_nu_nophase.insert({"IsNu", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::is_not_nu, &CUT, &TCUT)});
    }
    /// END GUNDAM VARIABLES ///////////////////////////////////////////////////////////////////////////////////////////////////
    //vars_selected_nu_nophase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::category, &CUT, &TCUT)});
    //vars_selected_nu_nophase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_muon_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_muon_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_leading_photon_cosphi", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_leading_photon_cosphi", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_subleading_photon_cosphi", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_subleading_photon_cosphi", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_pi0_mass", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_visible_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::visible_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_visible_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::visible_energy, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_nu_nophase.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_NoPhaseCuts", vars_selected_nu_nophase, false);
    */

    /*
    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_nophase;
    vars_selected_cos_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    /// GUNDAM VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////////////
    vars_selected_cos_nophase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::cut_type, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"IsNu", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::is_not_nu, &CUT, &TCUT)});
    if(std::string(argv[4]) == "sim")
      vars_selected_cos_nophase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::is_not_data, &CUT, &TCUT)});
    else
      vars_selected_cos_nophase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::is_data, &CUT, &TCUT)});
    /// END GUNDAM VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////////
    //vars_selected_cos_nophase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::category, &CUT, &TCUT)});
    //vars_selected_cos_nophase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_muon_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_muon_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_leading_photon_cosphi", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_leading_photon_cosphi", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_subleading_photon_cosphi", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_subleading_photon_cosphi", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_cosphi, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_pi0_mass", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_visible_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_nophase::visible_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_visible_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::visible_energy, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_cos_nophase.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_NoPhaseCuts", vars_selected_cos_nophase, false);
    */

    /*
    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_nophase;
    vars_purity_nophase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_nophase::category, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"one_muon_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::one_muon_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"zero_charged_pions_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::zero_charged_pions_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"two_or_three_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::two_or_three_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::all_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_nophase.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    if constexpr(WRITE_PURITY_TREES)
		  analysis.AddTree("Purity_NoPhaseCuts", vars_purity_nophase, false);
    */

    /*
    #undef SIGCUT
    #define SIGCUT cuts::ccpi0ana_nophase::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_nophase;
    vars_signal_nophase.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    /// GUNDAM VARIABLES /////////////////////////////////////////////////////////////////////////////////////////////////
    vars_signal_nophase.insert({"CutType", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::cut_type, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"IsData", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::is_not_data, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"IsNu", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::is_nu, &SIGCUT, &SIGCUT)});
    /// END GUNDAM VARIABLES /////////////////////////////////////////////////////////////////////////////////////////////
    //vars_signal_nophase.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    //vars_signal_nophase.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    //vars_signal_nophase.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    //vars_signal_nophase.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    //vars_signal_nophase.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    //vars_signal_nophase.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::category, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_leading_photon_cosphi", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_leading_photon_cosphi, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_subleading_photon_cosphi", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_subleading_photon_cosphi, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_nophase::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"one_muon_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::one_muon_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"zero_charged_pions_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::zero_charged_pions_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"two_or_three_photons_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::two_or_three_photons_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_vertex_x", SpineVar<TTYPE,TTYPE>(&vars::vertex_x, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_vertex_y", SpineVar<TTYPE,TTYPE>(&vars::vertex_y, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"true_vertex_z", SpineVar<TTYPE,TTYPE>(&vars::vertex_z, &SIGCUT, &SIGCUT)});
    vars_signal_nophase.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_nophase::all_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_NoPhaseCuts", vars_signal_nophase, true);
    */
    
    /**
     * @brief Add variabls for selected interactions (traditional) to the analysis.
     * @details This adds a set of variables to the analysis by creating a map
     * of variable names and SpillMultiVars that provide the functionality to
     * create the variables.  These names are used in the TTree that is created
     * by the Tree class to store the results of the analysis.
     */

    /*
    #undef CUT
    #define CUT cuts::ccpi0ana_trad::all_1mu0pi2gamma_cut
    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu_trad;
    vars_selected_nu_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_trad.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::is_signal_mc, &CUT, &TCUT)}); // GUNDA
    vars_selected_nu_trad.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_nu_trad.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::category, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::category_topology, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_muon_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_muon_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_pi0_mass", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_mass, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_visible_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::visible_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_visible_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::visible_energy, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_nu_trad.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedNu_TradCuts", vars_selected_nu_trad, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos_trad;
    vars_selected_cos_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::cut_type, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::is_signal_mc, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_selected_cos_trad.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::category, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::category_topology, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_muon_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_muon_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_muon_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_muon_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_pi0_leading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_pi0_leading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_pi0_subleading_photon_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_pi0_subleading_photon_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_conv_dist, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_pi0_momentum_mag", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_momentum_mag, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_pi0_beam_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_beam_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_pi0_photons_costheta", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_photons_costheta, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_pi0_mass", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_pi0_mass", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::pi0_mass, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_visible_energy", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_trad::visible_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_visible_energy", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::visible_energy, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_cos_trad.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedCos_TradCuts", vars_selected_cos_trad, false);


    #undef TCUT
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_trad;
    vars_purity_trad.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_trad::category, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"track_containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::track_containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"one_muon_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::one_muon_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"zero_charged_pions_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::zero_charged_pions_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"two_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::two_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"two_or_three_photons_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::two_or_three_photons_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"topology_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::topological_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"pi0_mass_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::pi0_mass_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"all_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::all_1mu0pi2gamma_cut), &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_trad.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    if constexpr(WRITE_PURITY_TREES)
		  analysis.AddTree("Purity_TradCuts", vars_purity_trad, false);

    #undef SIGCUT
    #define SIGCUT cuts::ccpi0ana_trad::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_signal_trad;
    vars_signal_trad.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"CutType", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::cut_type, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_trad.insert({"IsSignal", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::is_signal_mc, &SIGCUT, &SIGCUT)}); // GUNDAM
    vars_signal_trad.insert({"IsData", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::is_not_data, &SIGCUT, &SIGCUT)}); // GUNDAM 
    vars_signal_trad.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::category, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"category_topology", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::category_topology, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_muon_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::muon_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_muon_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::muon_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_pi0_leading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_pi0_leading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::pi0_leading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_pi0_subleading_photon_energy", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_energy, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_pi0_subleading_photon_conv_dist", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::pi0_subleading_photon_conv_dist, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_pi0_momentum_mag", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::pi0_momentum_mag, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_pi0_beam_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::pi0_beam_costheta, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_pi0_photons_costheta", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_pi0_mass", SpineVar<TTYPE,TTYPE>(&vars::ccpi0ana_trad::pi0_mass, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"track_containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::track_containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"one_muon_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::one_muon_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"zero_charged_pions_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::zero_charged_pions_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"two_or_three_photons_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::two_or_three_photons_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"topology_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::topological_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"pi0_mass_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::pi0_mass_cut), &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_vertex_x", SpineVar<TTYPE,TTYPE>(&vars::vertex_x, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_vertex_y", SpineVar<TTYPE,TTYPE>(&vars::vertex_y, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"true_vertex_z", SpineVar<TTYPE,TTYPE>(&vars::vertex_z, &SIGCUT, &SIGCUT)});
    vars_signal_trad.insert({"all_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::ccpi0ana_trad::all_1mu0pi2gamma_cut), &SIGCUT, &SIGCUT)});
    analysis.AddTree("Signal_TradCuts", vars_signal_trad, true);
    */

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
