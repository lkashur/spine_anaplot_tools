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
#include "include/srvar.h"
#include "include/analysis.h"
//#include "include/analysis_gundam.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void pi02024_gundam()
{
    ana::Analysis analysis("pi02024_output");

    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/fall2024/collonly_v2b/flat/*.root");
    //ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/fall2024/nominal/flat/*.root");
    analysis.AddLoader("mc", &mc, true);

    /**
     * @brief Add a set of variables for selected interactions relevant to the GUNDAM cross section analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::pi02024::all_1mu0pi2gamma_cut
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu;
    vars_selected_nu.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::pi02024::cut_type, &CUT, &TCUT)});
    vars_selected_nu.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::pi02024::is_signal, &CUT, &TCUT)});
    vars_selected_nu.insert({"IsData", SrVar<RTYPE,RTYPE>(&vars::pi02024::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu.insert({"RecoMuonP", SpineVar<RTYPE,RTYPE>(&vars::pi02024::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu.insert({"TrueMuonP", SpineVar<TTYPE,RTYPE>(&vars::pi02024::muon_momentum_mag, &CUT, &TCUT)});
    vars_selected_nu.insert({"RecoMuonCos", SpineVar<RTYPE,RTYPE>(&vars::pi02024::muon_beam_costheta, &CUT, &TCUT)});
    vars_selected_nu.insert({"TrueMuonCos", SpineVar<TTYPE,RTYPE>(&vars::pi02024::muon_beam_costheta, &CUT, &TCUT)});
    
    analysis.AddTree("selectedNu", vars_selected_nu, false);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
