/**
 * @file example_sbnd.C
 * @brief A basic example analysis macro for demonstrating the use of the
 * spine_anaplot_tools/cafana package on SBND simulation
 * @details This macro demonstrates how to use the spine_anaplot_tools/cafana
 * package to perform a basic analysis of a CAFAna file. The macro configures
 * some basic variables and cuts, then runs the analysis over a single sample.
 * @author mueller@fnal.gov
*/
#include "include/variables.h"
#include "include/muon2024/variables_muon2024.h"
#include "include/cuts.h"
#include "include/muon2024/cuts_muon2024.h"
#include "include/preprocessor.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void example_sbnd()
{
    /**
     * @brief Create an instance of the Analysis class.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefore the name of the output file, is specified as an argument
     * to the constructor.
     */
    ana::Analysis analysis("nuexample_sbnd_1muNp");

    /**
     * @brief Add a sample to the analysis.
     * @details This adds a sample to the analysis by creating a SpectrumLoader
     * object and adding it to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     */
    ana::SpectrumLoader sbnd("/pnfs/icarus/persistent/users/mueller/sbnd/larcv_sbnd_bnb_cosmics_spine.flat.root");
    analysis.AddLoader("sbnd", &sbnd, true);
    
    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::muon2024::all_1muNp_cut
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu;
    vars_selected_nu.insert({"nu_id", ana::SpillMultiVar(SPINEVAR_RT(vars::neutrino_id, CUT, TCUT))});
    vars_selected_nu.insert({"baseline", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_baseline, CUT, TCUT))});
    vars_selected_nu.insert({"pdg", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_pdg, CUT, TCUT))});
    vars_selected_nu.insert({"cc", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_cc, CUT, TCUT))});
    vars_selected_nu.insert({"category", ana::SpillMultiVar(SPINEVAR_RT(vars::muon2024::category, CUT, TCUT))});
    vars_selected_nu.insert({"interaction_mode", ana::SpillMultiVar(SPINEVAR_RT(vars::neutrino_interaction_mode, CUT, TCUT))});
    vars_selected_nu.insert({"true_edep", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_energy, CUT, TCUT))});
    vars_selected_nu.insert({"reco_edep", ana::SpillMultiVar(SPINEVAR_RR(vars::visible_energy, CUT, TCUT))});
    vars_selected_nu.insert({"true_muon_x", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_x, CUT, TCUT))});
    vars_selected_nu.insert({"reco_muon_x", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_x, CUT, TCUT))});
    vars_selected_nu.insert({"true_muon_y", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_y, CUT, TCUT))});
    vars_selected_nu.insert({"reco_muon_y", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_y, CUT, TCUT))});
    vars_selected_nu.insert({"true_muon_z", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_z, CUT, TCUT))});
    vars_selected_nu.insert({"reco_muon_z", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_z, CUT, TCUT))});
    vars_selected_nu.insert({"true_proton_x", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_x, CUT, TCUT))});
    vars_selected_nu.insert({"reco_proton_x", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_x, CUT, TCUT))});
    vars_selected_nu.insert({"true_proton_y", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_y, CUT, TCUT))});
    vars_selected_nu.insert({"reco_proton_y", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_y, CUT, TCUT))});
    vars_selected_nu.insert({"true_proton_z", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_z, CUT, TCUT))});
    vars_selected_nu.insert({"reco_proton_z", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_z, CUT, TCUT))});
    vars_selected_nu.insert({"true_tmuon", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_ke, CUT, TCUT))});
    vars_selected_nu.insert({"reco_tmuon", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_ke, CUT, TCUT))});
    vars_selected_nu.insert({"true_tproton", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_ke, CUT, TCUT))});
    vars_selected_nu.insert({"reco_tproton", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_ke, CUT, TCUT))});
    vars_selected_nu.insert({"true_ptmuon", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_pt, CUT, TCUT))});
    vars_selected_nu.insert({"reco_ptmuon", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_pt, CUT, TCUT))});
    vars_selected_nu.insert({"true_ptproton", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_pt, CUT, TCUT))});
    vars_selected_nu.insert({"reco_ptproton", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_pt, CUT, TCUT))});
    vars_selected_nu.insert({"true_theta_mu", ana::SpillMultiVar(SPINEVAR_RT(vars::muon_polar_angle, CUT, TCUT))});
    vars_selected_nu.insert({"reco_theta_mu", ana::SpillMultiVar(SPINEVAR_RR(vars::muon_polar_angle, CUT, TCUT))});
    vars_selected_nu.insert({"true_phi_mu", ana::SpillMultiVar(SPINEVAR_RT(vars::muon_azimuthal_angle, CUT, TCUT))});
    vars_selected_nu.insert({"reco_phi_mu", ana::SpillMultiVar(SPINEVAR_RR(vars::muon_azimuthal_angle, CUT, TCUT))});
    vars_selected_nu.insert({"true_opening_angle", ana::SpillMultiVar(SPINEVAR_RT(vars::muon2024::opening_angle, CUT, TCUT))});
    vars_selected_nu.insert({"reco_opening_angle", ana::SpillMultiVar(SPINEVAR_RR(vars::muon2024::opening_angle, CUT, TCUT))});
    vars_selected_nu.insert({"true_dpT", ana::SpillMultiVar(SPINEVAR_RT(vars::interaction_pt, CUT, TCUT))});
    vars_selected_nu.insert({"reco_dpT", ana::SpillMultiVar(SPINEVAR_RR(vars::interaction_pt, CUT, TCUT))});
    vars_selected_nu.insert({"true_dphiT", ana::SpillMultiVar(SPINEVAR_RT(vars::phiT, CUT, TCUT))});
    vars_selected_nu.insert({"reco_dphiT", ana::SpillMultiVar(SPINEVAR_RR(vars::phiT, CUT, TCUT))});
    vars_selected_nu.insert({"true_edalphaT", ana::SpillMultiVar(SPINEVAR_RT(vars::alphaT, CUT, TCUT))});
    vars_selected_nu.insert({"reco_edalphaT", ana::SpillMultiVar(SPINEVAR_RR(vars::alphaT, CUT, TCUT))});
    vars_selected_nu.insert({"true_vertex_x", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_x, CUT, TCUT))});
    vars_selected_nu.insert({"reco_vertex_x", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_x, CUT, TCUT))});
    vars_selected_nu.insert({"true_vertex_y", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_y, CUT, TCUT))});
    vars_selected_nu.insert({"reco_vertex_y", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_y, CUT, TCUT))});
    vars_selected_nu.insert({"true_vertex_z", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_z, CUT, TCUT))});
    vars_selected_nu.insert({"reco_vertex_z", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_z, CUT, TCUT))});
    vars_selected_nu.insert({"muon_muon_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_muon_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"muon_pion_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_pion_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"muon_proton_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_proton_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"muon_mip_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_mip_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"muon_hadron_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_hadron_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"proton_muon_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_muon_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"proton_pion_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_pion_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"proton_proton_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_proton_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"proton_mip_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_mip_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"proton_hadron_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_hadron_softmax, CUT, TCUT))});
    vars_selected_nu.insert({"flash_time", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_time, CUT, TCUT))});
    vars_selected_nu.insert({"flash_total", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_total_pe, CUT, TCUT))});
    vars_selected_nu.insert({"flash_hypothesis", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_hypothesis, CUT, TCUT))});

    analysis.AddTree("selectedNu", vars_selected_nu, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos;
    vars_selected_cos.insert({"nu_id", ana::SpillMultiVar(SPINEVAR_RT(vars::neutrino_id, CUT, TCUT))});
    vars_selected_cos.insert({"baseline", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_baseline, CUT, TCUT))});
    vars_selected_cos.insert({"pdg", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_pdg, CUT, TCUT))});
    vars_selected_cos.insert({"cc", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_cc, CUT, TCUT))});
    vars_selected_cos.insert({"category", ana::SpillMultiVar(SPINEVAR_RT(vars::muon2024::category, CUT, TCUT))});
    vars_selected_cos.insert({"interaction_mode", ana::SpillMultiVar(SPINEVAR_RT(vars::neutrino_interaction_mode, CUT, TCUT))});
    vars_selected_cos.insert({"true_edep", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_energy, CUT, TCUT))});
    vars_selected_cos.insert({"reco_edep", ana::SpillMultiVar(SPINEVAR_RR(vars::visible_energy, CUT, TCUT))});
    vars_selected_cos.insert({"true_muon_x", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_x, CUT, TCUT))});
    vars_selected_cos.insert({"reco_muon_x", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_x, CUT, TCUT))});
    vars_selected_cos.insert({"true_muon_y", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_y, CUT, TCUT))});
    vars_selected_cos.insert({"reco_muon_y", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_y, CUT, TCUT))});
    vars_selected_cos.insert({"true_muon_z", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_z, CUT, TCUT))});
    vars_selected_cos.insert({"reco_muon_z", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_z, CUT, TCUT))});
    vars_selected_cos.insert({"true_proton_x", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_x, CUT, TCUT))});
    vars_selected_cos.insert({"reco_proton_x", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_x, CUT, TCUT))});
    vars_selected_cos.insert({"true_proton_y", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_y, CUT, TCUT))});
    vars_selected_cos.insert({"reco_proton_y", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_y, CUT, TCUT))});
    vars_selected_cos.insert({"true_proton_z", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_z, CUT, TCUT))});
    vars_selected_cos.insert({"reco_proton_z", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_z, CUT, TCUT))});
    vars_selected_cos.insert({"true_tmuon", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_ke, CUT, TCUT))});
    vars_selected_cos.insert({"reco_tmuon", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_ke, CUT, TCUT))});
    vars_selected_cos.insert({"true_tproton", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_ke, CUT, TCUT))});
    vars_selected_cos.insert({"reco_tproton", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_ke, CUT, TCUT))});
    vars_selected_cos.insert({"true_ptmuon", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_pt, CUT, TCUT))});
    vars_selected_cos.insert({"reco_ptmuon", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_pt, CUT, TCUT))});
    vars_selected_cos.insert({"true_ptproton", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_pt, CUT, TCUT))});
    vars_selected_cos.insert({"reco_ptproton", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_pt, CUT, TCUT))});
    vars_selected_cos.insert({"true_theta_mu", ana::SpillMultiVar(SPINEVAR_RT(vars::muon_polar_angle, CUT, TCUT))});
    vars_selected_cos.insert({"reco_theta_mu", ana::SpillMultiVar(SPINEVAR_RR(vars::muon_polar_angle, CUT, TCUT))});
    vars_selected_cos.insert({"true_phi_mu", ana::SpillMultiVar(SPINEVAR_RT(vars::muon_azimuthal_angle, CUT, TCUT))});
    vars_selected_cos.insert({"reco_phi_mu", ana::SpillMultiVar(SPINEVAR_RR(vars::muon_azimuthal_angle, CUT, TCUT))});
    vars_selected_cos.insert({"true_opening_angle", ana::SpillMultiVar(SPINEVAR_RT(vars::muon2024::opening_angle, CUT, TCUT))});
    vars_selected_cos.insert({"reco_opening_angle", ana::SpillMultiVar(SPINEVAR_RR(vars::muon2024::opening_angle, CUT, TCUT))});
    vars_selected_cos.insert({"true_dpT", ana::SpillMultiVar(SPINEVAR_RT(vars::interaction_pt, CUT, TCUT))});
    vars_selected_cos.insert({"reco_dpT", ana::SpillMultiVar(SPINEVAR_RR(vars::interaction_pt, CUT, TCUT))});
    vars_selected_cos.insert({"true_dphiT", ana::SpillMultiVar(SPINEVAR_RT(vars::phiT, CUT, TCUT))});
    vars_selected_cos.insert({"reco_dphiT", ana::SpillMultiVar(SPINEVAR_RR(vars::phiT, CUT, TCUT))});
    vars_selected_cos.insert({"true_edalphaT", ana::SpillMultiVar(SPINEVAR_RT(vars::alphaT, CUT, TCUT))});
    vars_selected_cos.insert({"reco_edalphaT", ana::SpillMultiVar(SPINEVAR_RR(vars::alphaT, CUT, TCUT))});
    vars_selected_cos.insert({"true_vertex_x", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_x, CUT, TCUT))});
    vars_selected_cos.insert({"reco_vertex_x", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_x, CUT, TCUT))});
    vars_selected_cos.insert({"true_vertex_y", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_y, CUT, TCUT))});
    vars_selected_cos.insert({"reco_vertex_y", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_y, CUT, TCUT))});
    vars_selected_cos.insert({"true_vertex_z", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_z, CUT, TCUT))});
    vars_selected_cos.insert({"reco_vertex_z", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_z, CUT, TCUT))});
    vars_selected_cos.insert({"muon_muon_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_muon_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"muon_pion_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_pion_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"muon_proton_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_proton_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"muon_mip_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_mip_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"muon_hadron_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_hadron_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"proton_muon_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_muon_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"proton_pion_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_pion_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"proton_proton_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_proton_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"proton_mip_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_mip_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"proton_hadron_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_hadron_softmax, CUT, TCUT))});
    vars_selected_cos.insert({"flash_time", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_time, CUT, TCUT))});
    vars_selected_cos.insert({"flash_total", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_total_pe, CUT, TCUT))});
    vars_selected_cos.insert({"flash_hypothesis", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_hypothesis, CUT, TCUT))});

    analysis.AddTree("selectedCos", vars_selected_cos, false);

    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_purity_nu;
    vars_purity_nu.insert({"nu_id", ana::SpillMultiVar(SPINEVAR_RT(vars::neutrino_id, cuts::no_cut, TCUT))});
    vars_purity_nu.insert({"category", ana::SpillMultiVar(SPINEVAR_RT(vars::muon2024::category, cuts::no_cut, TCUT))});
    vars_purity_nu.insert({"fiducial_cut", ana::SpillMultiVar(SPINEVAR_RR(cuts::fiducial_cut, cuts::no_cut, TCUT))});
    vars_purity_nu.insert({"containment_cut", ana::SpillMultiVar(SPINEVAR_RR(cuts::containment_cut, cuts::no_cut, TCUT))});
    vars_purity_nu.insert({"flash_cut_bnb", ana::SpillMultiVar(SPINEVAR_RR(cuts::flash_cut_bnb, cuts::no_cut, TCUT))});
    vars_purity_nu.insert({"has_no_charged_pions", ana::SpillMultiVar(SPINEVAR_RR(cuts::no_charged_pions, cuts::no_cut, TCUT))});
    vars_purity_nu.insert({"has_no_showers", ana::SpillMultiVar(SPINEVAR_RR(cuts::no_showers, cuts::no_cut, TCUT))});
    vars_purity_nu.insert({"has_single_muon", ana::SpillMultiVar(SPINEVAR_RR(cuts::has_single_muon, cuts::no_cut, TCUT))});
    vars_purity_nu.insert({"has_multiple_protons", ana::SpillMultiVar(SPINEVAR_RR(cuts::has_nonzero_protons, cuts::no_cut, TCUT))});

    analysis.AddTree("purityNu", vars_purity_nu, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_purity_cos;
    vars_purity_cos.insert({"nu_id", ana::SpillMultiVar(SPINEVAR_RT(vars::neutrino_id, cuts::no_cut, TCUT))});
    vars_purity_cos.insert({"category", ana::SpillMultiVar(SPINEVAR_RT(vars::muon2024::category, cuts::no_cut, TCUT))});
    vars_purity_cos.insert({"fiducial_cut", ana::SpillMultiVar(SPINEVAR_RR(cuts::fiducial_cut, cuts::no_cut, TCUT))});
    vars_purity_cos.insert({"containment_cut", ana::SpillMultiVar(SPINEVAR_RR(cuts::containment_cut, cuts::no_cut, TCUT))});
    vars_purity_cos.insert({"flash_cut_bnb", ana::SpillMultiVar(SPINEVAR_RR(cuts::flash_cut_bnb, cuts::no_cut, TCUT))});
    vars_purity_cos.insert({"has_no_charged_pions", ana::SpillMultiVar(SPINEVAR_RR(cuts::no_charged_pions, cuts::no_cut, TCUT))});
    vars_purity_cos.insert({"has_no_showers", ana::SpillMultiVar(SPINEVAR_RR(cuts::no_showers, cuts::no_cut, TCUT))});
    vars_purity_cos.insert({"has_single_muon", ana::SpillMultiVar(SPINEVAR_RR(cuts::has_single_muon, cuts::no_cut, TCUT))});
    vars_purity_cos.insert({"has_multiple_protons", ana::SpillMultiVar(SPINEVAR_RR(cuts::has_nonzero_protons, cuts::no_cut, TCUT))});

    analysis.AddTree("purityCos", vars_purity_cos, false);

    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define SIGCUT cuts::muon2024::signal_1muNp
    std::map<std::string, ana::SpillMultiVar> vars_signal;
    vars_signal.insert({"nu_id", ana::SpillMultiVar(SPINEVAR_TT(vars::neutrino_id, SIGCUT))});
    vars_signal.insert({"baseline", ana::SpillMultiVar(SPINEVAR_TT(vars::true_neutrino_baseline, SIGCUT))});
    vars_signal.insert({"pdg", ana::SpillMultiVar(SPINEVAR_TT(vars::true_neutrino_pdg, SIGCUT))});
    vars_signal.insert({"cc", ana::SpillMultiVar(SPINEVAR_TT(vars::true_neutrino_cc, SIGCUT))});
    vars_signal.insert({"category", ana::SpillMultiVar(SPINEVAR_TT(vars::muon2024::category, SIGCUT))});
    vars_signal.insert({"interaction_mode", ana::SpillMultiVar(SPINEVAR_TT(vars::neutrino_interaction_mode, SIGCUT))});
    vars_signal.insert({"true_edep", ana::SpillMultiVar(SPINEVAR_TT(vars::true_neutrino_energy, SIGCUT))});
    vars_signal.insert({"true_muon_x", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_end_x, SIGCUT))});
    vars_signal.insert({"true_muon_y", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_end_y, SIGCUT))});
    vars_signal.insert({"true_muon_z", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_end_z, SIGCUT))});
    vars_signal.insert({"true_proton_x", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_end_x, SIGCUT))});
    vars_signal.insert({"true_proton_y", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_end_y, SIGCUT))});
    vars_signal.insert({"true_proton_z", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_end_z, SIGCUT))});
    vars_signal.insert({"true_tmuon", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_ke, SIGCUT))});
    vars_signal.insert({"true_tproton", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_ke, SIGCUT))});
    vars_signal.insert({"true_ptmuon", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_pt, SIGCUT))});
    vars_signal.insert({"true_ptproton", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_pt, SIGCUT))});
    vars_signal.insert({"true_theta_mu", ana::SpillMultiVar(SPINEVAR_TT(vars::muon_polar_angle, SIGCUT))});
    vars_signal.insert({"true_phi_mu", ana::SpillMultiVar(SPINEVAR_TT(vars::muon_azimuthal_angle, SIGCUT))});
    vars_signal.insert({"true_opening_angle", ana::SpillMultiVar(SPINEVAR_TT(vars::muon2024::opening_angle, SIGCUT))});
    vars_signal.insert({"true_dpT", ana::SpillMultiVar(SPINEVAR_TT(vars::interaction_pt, SIGCUT))});
    vars_signal.insert({"true_dphiT", ana::SpillMultiVar(SPINEVAR_TT(vars::phiT, SIGCUT))});
    vars_signal.insert({"true_edalphaT", ana::SpillMultiVar(SPINEVAR_TT(vars::alphaT, SIGCUT))});
    vars_signal.insert({"true_vertex_x", ana::SpillMultiVar(SPINEVAR_TT(vars::vertex_x, SIGCUT))});
    vars_signal.insert({"true_vertex_y", ana::SpillMultiVar(SPINEVAR_TT(vars::vertex_y, SIGCUT))});
    vars_signal.insert({"true_vertex_z", ana::SpillMultiVar(SPINEVAR_TT(vars::vertex_z, SIGCUT))});
    vars_signal.insert({"muon_muon_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_muon_muon_softmax, SIGCUT))});
    vars_signal.insert({"muon_pion_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_muon_pion_softmax, SIGCUT))});
    vars_signal.insert({"muon_proton_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_muon_proton_softmax, SIGCUT))});
    vars_signal.insert({"muon_mip_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_muon_mip_softmax, SIGCUT))});
    vars_signal.insert({"muon_hadron_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_muon_hadron_softmax, SIGCUT))});
    vars_signal.insert({"proton_muon_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_proton_muon_softmax, SIGCUT))});
    vars_signal.insert({"proton_pion_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_proton_pion_softmax, SIGCUT))});
    vars_signal.insert({"proton_proton_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_proton_proton_softmax, SIGCUT))});
    vars_signal.insert({"proton_mip_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_proton_mip_softmax, SIGCUT))});
    vars_signal.insert({"proton_hadron_softmax", ana::SpillMultiVar(SPINEVAR_TR(vars::leading_proton_hadron_softmax, SIGCUT))});
    vars_signal.insert({"fiducial_cut", ana::SpillMultiVar(SPINEVAR_TR(cuts::fiducial_cut, SIGCUT))});
    vars_signal.insert({"containment_cut", ana::SpillMultiVar(SPINEVAR_TR(cuts::containment_cut, SIGCUT))});
    vars_signal.insert({"flash_cut_bnb", ana::SpillMultiVar(SPINEVAR_TR(cuts::flash_cut_bnb, SIGCUT))});
    vars_signal.insert({"no_charged_pions", ana::SpillMultiVar(SPINEVAR_TR(cuts::no_charged_pions, SIGCUT))});
    vars_signal.insert({"no_showers", ana::SpillMultiVar(SPINEVAR_TR(cuts::no_showers, SIGCUT))});
    vars_signal.insert({"has_single_muon", ana::SpillMultiVar(SPINEVAR_TR(cuts::has_single_muon, SIGCUT))});
    vars_signal.insert({"has_multiple_protons", ana::SpillMultiVar(SPINEVAR_TR(cuts::has_nonzero_protons, SIGCUT))});
    
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