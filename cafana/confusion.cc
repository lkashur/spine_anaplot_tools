/**
 * @file confusion.C
 * @brief The main analysis macro for the pi0ana analysis on ICARUS Monte
 * Carlo simulation.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author lkashur@colostate.edu
*/
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PRIMARYFUNC pvars::primary_classification
#define PIDFUNC pvars::pid
//#define PIDFUNC pvars::custom_pid                                                                                                                                                    
#define PROTON_BINDING_ENERGY 30.9 // MeV                                                                                                                                              
#define BEAM_IS_NUMI false
#define WRITE_PURITY_TREES false

#include "include/mctruth.h"
#include "include/variables.h"
#include "include/cuts.h"

#include "include/pi0ana/variables_ccpi0ana_phase.h"
#include "include/pi0ana/cuts_ccpi0ana_phase.h"

#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

int main()
{

    // Configure analysis
    ana::Analysis analysis("icarus_bnb_ccpi0_confusion");
    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/mixed/simulation/nominal/input*flat.root");
    analysis.AddLoader("mc", &mc, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define SIGCUT cuts::ccpi0ana_phase::signal_1mu0pi1pi0
    std::map<std::string, ana::SpillMultiVar> vars_eff_pid_confusion_phase;
    auto true_primary_particle = [](const TTYPEP & p) { return p.is_primary && PIDFUNC(p) >=0; };
    auto true_pid = [](const TTYPEP & p) -> double { return PIDFUNC(p); };
    auto reco_pid = [](const RTYPEP & p) -> double { return PIDFUNC(p); };
    vars_eff_pid_confusion_phase.insert({"true_pid", SpineVar<TTYPEP,TTYPEP,TTYPE>(true_pid, true_primary_particle, &SIGCUT)});
    vars_eff_pid_confusion_phase.insert({"reco_pid", SpineVar<RTYPEP,TTYPEP,TTYPE>(reco_pid, true_primary_particle, &SIGCUT)});
    vars_eff_pid_confusion_phase.insert({"reco_primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), true_primary_particle, &SIGCUT)});
    analysis.AddTree("EffPIDConfusion_PhaseCuts", vars_eff_pid_confusion_phase, true);

    /*
    std::map<std::string, ana::SpillMultiVar> vars_eff_primary_confusion_phase;
    auto true_particle = [](const TTYPEP & p) { return p.pid >=0; };
    vars_eff_primary_confusion_phase.insert({"true_pid", SpineVar<TTYPEP,TTYPEP,TTYPE>(true_pid, true_particle, &SIGCUT)});
    vars_eff_primary_confusion_phase.insert({"reco_pid", SpineVar<RTYPEP,TTYPEP,TTYPE>(reco_pid, true_particle, &SIGCUT)});
    vars_eff_primary_confusion_phase.insert({"true_primary", SpineVar<TTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), true_particle, &SIGCUT)});
    vars_eff_primary_confusion_phase.insert({"reco_primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), true_particle, &SIGCUT)});
    analysis.AddTree("EffPrimaryConfusion_PhaseCuts", vars_eff_primary_confusion_phase, true);
    */
    
    #define CUT cuts::ccpi0ana_phase::all_cut
    std::map<std::string, ana::SpillMultiVar> vars_pur_pid_confusion_phase;
    auto reco_primary_particle = [](const RTYPEP & p) { return p.is_primary && PIDFUNC(p) >=0; };
    vars_pur_pid_confusion_phase.insert({"reco_pid", SpineVar<RTYPEP,RTYPEP,RTYPE>(reco_pid, reco_primary_particle, &CUT)});
    //vars_pur_pid_confusion_phase.insert({"true_pid", SpineVar<TTYPEP,RTYPEP,RTYPE>(true_pid, reco_primary_particle, &CUT)});
    //vars_pur_pid_confusion_phase.insert({"true_primary", SpineVar<TTYPEP,RTYPEP,RTYPE>(WRAP_BOOL(pcuts::is_primary), reco_primary_particle, &CUT)});
    analysis.AddTree("PurPIDConfusion_PhaseCuts", vars_pur_pid_confusion_phase, true);
    

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
