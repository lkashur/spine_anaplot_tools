/**
 * @file energycalib.C
 * @brief The main analysis macro for the muon2024 analysis on ICARUS Monte
 * Carlo simulation.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author mueller@fnal.gov
*/
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PIDFUNC pvars::pid
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false

#include "include/mctruth.h"
#include "include/variables.h"
#include "include/muon2024/variables_muon2024.h"
#include "include/cuts.h"
#include "include/muon2024/cuts_muon2024.h"
#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void energycalib()
{

    // Configure analysis 
    ana::Analysis analysis("energycalib_bnb_cv_05_march_2025");
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/lkashur/bnb_nu_cosmic_cv_merged/flatcaf*.root"); // old CV
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/production/simulation/nominal/flat/input*.flat.root"); // updated CV
    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/test/nominal_v4/flat/input*.flat.root"); // BNB nu + cosmic CV, upstream calibration fixed
    analysis.AddLoader("mc", &mc, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    std::map<std::string, ana::SpillMultiVar> vars_muon;
    auto primary_muon = [](const TTYPEP & p) { return p.pid == 2 && p.is_primary && (0.001*p.t < 1.6 && 0.001*p.t > 0);};
    auto primary_muon_reco = [](const RTYPEP & p) -> double { return p.pid == 2 && p.is_primary; };
    vars_muon.insert({"true_muon_ke", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::ke, primary_muon, &cuts::no_cut)});
    vars_muon.insert({"reco_csda_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::csda_ke, primary_muon, &cuts::no_cut)});
    vars_muon.insert({"reco_calo_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::calo_ke, primary_muon, &cuts::no_cut)});
    vars_muon.insert({"reco_mcs_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::mcs_ke, primary_muon, &cuts::no_cut)});
    vars_muon.insert({"reco_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::ke, primary_muon, &cuts::no_cut)});
    vars_muon.insert({"true_is_contained", SpineVar<TTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_contained), primary_muon, &cuts::no_cut)});
    vars_muon.insert({"reco_is_contained", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_contained), primary_muon, &cuts::no_cut)});
    vars_muon.insert({"reco_is_primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), primary_muon, &cuts::no_cut)});
    vars_muon.insert({"reco_is_muon", SpineVar<RTYPEP,TTYPEP,TTYPE>(primary_muon_reco, primary_muon, &cuts::no_cut)});
    analysis.AddTree("Muon", vars_muon, true);

    std::map<std::string, ana::SpillMultiVar> vars_photon;
    auto primary_photon = [](const TTYPEP & p) { return p.pid == 0 && p.is_primary && p.is_contained && (0.001*p.t < 1.6 && 0.001*p.t > 0);};
    auto primary_photon_reco = [](const RTYPEP & p) -> double { return p.pid == 0 && p.is_primary; };
    vars_photon.insert({"true_photon_ke", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::ke, primary_photon, &cuts::no_cut)});
    vars_photon.insert({"reco_calo_ke_pre_corr", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::calo_ke_pre_corr, primary_photon, &cuts::no_cut)});
    vars_photon.insert({"reco_calo_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::calo_ke, primary_photon, &cuts::no_cut)});
    vars_photon.insert({"reco_is_contained", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_contained), primary_photon, &cuts::no_cut)});
    vars_photon.insert({"reco_is_primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), primary_photon, &cuts::no_cut)});
    vars_photon.insert({"reco_is_photon", SpineVar<RTYPEP,TTYPEP,TTYPE>(primary_photon_reco, primary_photon, &cuts::no_cut)});
    analysis.AddTree("Photon", vars_photon, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
