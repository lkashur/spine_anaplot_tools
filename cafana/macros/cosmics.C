/**
 * @file cosmics.C
 * @brief An analysis macro for basic cosmic muon studies.
 * @details This macro is used to perform a basic analysis of cosmic muons in
 * the ICARUS TPCs. The macro configures some basic variables and cuts, then
 * runs the analysis over a single sample of off-beam data. The results of the
 * analysis are stored in a ROOT file that can be used to create plots and
 * perform further studies.
 * @author mueller@fnal.gov
*/
#include "include/variables.h"
#include "include/cuts.h"
#include "include/cosmics/cuts_cosmics.h"
#include "include/preprocessor.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

void cosmics()
{
    /**
     * @brief Create an instance of the Analysis class.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefore the name of the output file, is specified as an argument
     * to the constructor.
     */
    ana::Analysis analysis("cosmics");

    /**
     * @brief Add samples to the analysis.
     * @details This adds samples to the analysis by creating SpectrumLoader
     * objects and adding them to the Analysis class. Here, we add an on-beam
     * MC sample and an off-beam data sample to the analysis. These have the
     * same triggering logic, so we can make meaningful comparisons between
     * the two samples.
     */
    ana::SpectrumLoader onbeam("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/mc_v09_84_00_01/flat/*.root");
    analysis.AddLoader("mc_onbeam", &onbeam, true);

    ana::SpectrumLoader offbeam("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/data_v09_84_00_01/offbeam/flat/*.root");
    analysis.AddLoader("data_offbeam", &offbeam, false);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    std::map<std::string, ana::SpillMultiVar> vars;
    #define CUT cuts::cosmics::single_cosmic_muon_cut
    vars.insert({"muon_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_softmax, CUT, cuts::no_cut))});
    vars.insert({"proton_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_softmax, CUT, cuts::no_cut))});
    vars.insert({"mip_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_mip_softmax, CUT, cuts::no_cut))});
    vars.insert({"flash_time", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_time, CUT, cuts::no_cut))});
    vars.insert({"flash_total", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_total_pe, CUT, cuts::no_cut))});
    vars.insert({"flash_hypothesis", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_hypothesis, CUT, cuts::no_cut))});
    vars.insert({"true_containment", ana::SpillMultiVar(SPINEVAR_RT(vars::containment, CUT, cuts::no_cut))});
    vars.insert({"reco_containment", ana::SpillMultiVar(SPINEVAR_RR(vars::containment, CUT, cuts::no_cut))});
    vars.insert({"true_fiducial", ana::SpillMultiVar(SPINEVAR_RT(vars::fiducial, CUT, cuts::no_cut))});
    vars.insert({"reco_fiducial", ana::SpillMultiVar(SPINEVAR_RR(vars::fiducial, CUT, cuts::no_cut))});
    vars.insert({"true_muon_x", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_x, CUT, cuts::no_cut))});
    vars.insert({"reco_muon_x", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_x, CUT, cuts::no_cut))});
    vars.insert({"true_muon_y", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_y, CUT, cuts::no_cut))});
    vars.insert({"reco_muon_y", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_y, CUT, cuts::no_cut))});
    vars.insert({"true_muon_z", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_z, CUT, cuts::no_cut))});
    vars.insert({"reco_muon_z", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_z, CUT, cuts::no_cut))});
    vars.insert({"true_vertex_x", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_x, CUT, cuts::no_cut))});
    vars.insert({"reco_vertex_x", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_x, CUT, cuts::no_cut))});
    vars.insert({"true_vertex_y", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_y, CUT, cuts::no_cut))});
    vars.insert({"reco_vertex_y", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_y, CUT, cuts::no_cut))});
    vars.insert({"true_vertex_z", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_z, CUT, cuts::no_cut))});
    vars.insert({"reco_vertex_z", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_z, CUT, cuts::no_cut))});
    vars.insert({"true_tmuon", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_ke, CUT, cuts::no_cut))});
    vars.insert({"reco_tmuon", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_ke, CUT, cuts::no_cut))});
    vars.insert({"true_ptmuon", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_pt, CUT, cuts::no_cut))});
    vars.insert({"reco_ptmuon", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_pt, CUT, cuts::no_cut))});
    vars.insert({"true_theta_mu", ana::SpillMultiVar(SPINEVAR_RT(vars::muon_polar_angle, CUT, cuts::no_cut))});
    vars.insert({"reco_theta_mu", ana::SpillMultiVar(SPINEVAR_RR(vars::muon_polar_angle, CUT, cuts::no_cut))});
    vars.insert({"true_phi_mu", ana::SpillMultiVar(SPINEVAR_RT(vars::muon_azimuthal_angle, CUT, cuts::no_cut))});
    vars.insert({"reco_phi_mu", ana::SpillMultiVar(SPINEVAR_RR(vars::muon_azimuthal_angle, CUT, cuts::no_cut))});
    #undef CUT

    analysis.AddTree("cosmics", vars, false);

    /**
     * @brief Run the analysis on the specified samples.
     * @details This function runs the analysis on the samples that have been
     * added to the Analysis class. The results of the analysis are stored in
     * a ROOT file that can be used to create plots and perform further studies.
     */
    analysis.Go();
}