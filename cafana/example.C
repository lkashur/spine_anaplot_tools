/**
 * @file example.C
 * @brief A basic example analysis macro for demonstrating the use of the
 * spine_anaplot_tools/cafana package.
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

#include <algorithm>

namespace ana
{ 
    #define CUTNAME cuts::fiducial_containment_flash_cut_bnb
    VARDLP_RECO_TRUE(kTrueE, vars::true_neutrino_energy, CUTNAME);
    VARDLP_RECO_TRUE(kTrueL, vars::true_neutrino_baseline, CUTNAME);
    VARDLP_RECO_TRUE(kTruePDG, vars::true_neutrino_pdg, CUTNAME);
    VARDLP_RECO_TRUE(kTrueCC, vars::true_neutrino_cc, CUTNAME);
    VARDLP_RECO_RECO(kRecoE, vars::visible_energy, CUTNAME);
}

void example()
{
    /**
     * @brief Create an instance of the Analysis class.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefor the name of the output file, is specified as an argument to
     * the constructor.
     */
    ana::Analysis analysis("nuexample");

    /**
     * @brief Add a sample to the analysis.
     * @details This adds a sample to the analysis by creating a SpectrumLoader
     * object and adding it to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     */
    ana::SpectrumLoader var00("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var00_nominal.flat.root");
    analysis.AddLoader("nominal", &var00);
    
    /**
     * @brief Add a set of variables to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * vector of variable names and a vector with the corresponding
     * SpillMultiVars that produce the variables. These names are used in the
     * TTree that is created by the Tree class to store the results of the
     * analysis.
     */
    std::vector<std::string> names({"trueE", "trueL", "truePDG", "CC", "recoE"});
    std::vector<ana::SpillMultiVar> vars({ana::kTrueE, ana::kTrueL, ana::kTruePDG, ana::kTrueCC, ana::kRecoE});
    analysis.AddVars(names, vars);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}