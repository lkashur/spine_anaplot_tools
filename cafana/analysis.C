/**
 * @file profit.C
 * @brief ROOT macro to be used with CAFAna to run the selection and produce sBruce trees.
 * @author mueller@fnal.gov
*/
#include "include/variables.h"
#include "include/muon2024/variables_muon2024.h"
#include "include/cuts.h"
#include "include/muon2024/cuts_muon2024.h"
#include "include/preprocessor.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

#include <algorithm>

namespace ana
{ 
    VARDLP_RECO_TRUE(kTrueE, vars::true_neutrino_energy, cuts::muon2024::all_1mu1p_cut);
    VARDLP_RECO_TRUE(kTrueL, vars::true_neutrino_baseline, cuts::muon2024::all_1mu1p_cut);
    VARDLP_RECO_TRUE(kTruePDG, vars::true_neutrino_pdg, cuts::muon2024::all_1mu1p_cut);
    VARDLP_RECO_TRUE(kTrueCC, vars::true_neutrino_cc, cuts::muon2024::all_1mu1p_cut);
    VARDLP_RECO_RECO(kRecoE, vars::visible_energy, cuts::muon2024::all_1mu1p_cut);
    VARDLP_RECO_RECO(kMuonSoftmax, vars::leading_muon_softmax, cuts::muon2024::all_1mu1p_cut);
    VARDLP_RECO_RECO(kProtonSoftmax, vars::leading_proton_softmax, cuts::muon2024::all_1mu1p_cut);
}

void analysis()
{
    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var00_nominal.flat.root");    
    ana::Tree nutree("selectedNu", {"trueE", "trueL", "truePDG", "CC", "recoE", "muon_pid", "proton_pid"}, mc,
                    {ana::kTrueE, ana::kTrueL, ana::kTruePDG, ana::kTrueCC, ana::kRecoE, ana::kMuonSoftmax, ana::kProtonSoftmax},
                    ana::kNoSpillCut,true);
    mc.Go();
    TFile output("var00_nominal.root", "RECREATE");
    TDirectory * dir = output.mkdir("events");
    nutree.SaveTo(dir);
}
