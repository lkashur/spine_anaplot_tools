/**
 * @file profit.C
 * @brief ROOT macro to be used with CAFAna to run the selection and produce sBruce trees.
 * @author mueller@fnal.gov
*/
#include "include/variables.h"
#include "include/numu_variables.h"
#include "include/cuts.h"
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
    VARDLP_RECO_TRUE(kTrueE, vars::interaction::true_neutrino_energy, cuts::all_1mu1p_cut);
    VARDLP_RECO_TRUE(kTrueL, vars::interaction::true_neutrino_baseline, cuts::all_1mu1p_cut);
    VARDLP_RECO_TRUE(kTruePDG, vars::interaction::true_neutrino_pdg, cuts::all_1mu1p_cut);
    VARDLP_RECO_TRUE(kTrueCC, vars::interaction::true_neutrino_cc, cuts::all_1mu1p_cut);
    VARDLP_RECO_RECO(kRecoE, vars:interaction::visible_energy, cuts::all_1mu1p_cut);
    VARDLP_RECO_RECO(kMuonSoftmax, vars::interaction::muon_softmax, cuts::all_1mu1p_cut);
    VARDLP_RECO_RECO(kProtonSoftmax, vars::interaction::proton_softmax, cuts::all_1mu1p_cut);
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
