/**
 * @file analysis.C
 * @brief The main analysis macro for the muon2024 analysis.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
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
    VARDLP_RECO_TRUE(kTrueE, vars::true_neutrino_energy, cuts::muon2024::all_1muNp_cut);
    VARDLP_RECO_TRUE(kTrueL, vars::true_neutrino_baseline, cuts::muon2024::all_1muNp_cut);
    VARDLP_RECO_TRUE(kTruePDG, vars::true_neutrino_pdg, cuts::muon2024::all_1muNp_cut);
    VARDLP_RECO_TRUE(kTrueCC, vars::true_neutrino_cc, cuts::muon2024::all_1muNp_cut);
    VARDLP_RECO_RECO(kRecoE, vars::visible_energy, cuts::muon2024::all_1muNp_cut);
    VARDLP_RECO_RECO(kMuonSoftmax, vars::leading_muon_softmax, cuts::muon2024::all_1muNp_cut);
    VARDLP_RECO_RECO(kProtonSoftmax, vars::leading_proton_softmax, cuts::muon2024::all_1muNp_cut);
    VARDLP_RECO_RECO(kMIPSoftmax, vars::leading_muon_mip_softmax, cuts::muon2024::all_1muNp_cut);
}

void analysis()
{
    ana::Analysis analysis("muon2024_1muNp");

    ana::SpectrumLoader var00("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var00_nominal.flat.root");
    analysis.AddLoader("nominal", &var00);
    std::vector<std::string> names({"trueE", "trueL", "truePDG", "CC", "recoE", "muon_pid", "proton_pid", "mip_pid"});
    std::vector<ana::SpillMultiVar> vars({ana::kTrueE, ana::kTrueL, ana::kTruePDG, ana::kTrueCC, ana::kRecoE, ana::kMuonSoftmax, ana::kProtonSoftmax, ana::kMIPSoftmax});
    analysis.AddVars(names, vars);
    analysis.Go();
}
