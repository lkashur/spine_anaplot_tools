vector<string> GetGENIEMultisigmaKnobNames();

void spine2gundam()
{
  
  vector<string> genieMultisigmaNames = GetGENIEMultisigmaKnobNames();
 
  // Input
  string infile_string = "icarus_bnb_ccpi0_mc_onbeam_offbeam_vars_syst_delete.root";
  string selnutree_string = "/events/mc/SelectedNu_PhaseCuts";
  string selcostree_string = "/events/mc/SelectedCos_PhaseCuts";
  string potmc_string = "/events/mc/POT";
  string livetimemc_string = "events/mc/Livetime";
  string seloffbeamtree_string = "/events/offbeam/SelectedCos_PhaseCuts";
  string livetimeoffbeam_string = "events/offbeam/Livetime";
  string selonbeamtree_string = "/events/onbeam/SelectedNu_PhaseCuts";
  string potonbeam_string = "/events/onbeam/POT";
  string livetimeonbeam_string = "/events/onbeam/Livetime";
  string signaltree_string = "events/mc/Signal_PhaseCuts";
  
  // Output
  string outfile_mc_string = "icarus_bnb_ccpi0_mc_offbeam_syst_gundaminput.root"; // mc = CV + off-beam
  string outfile_data_string = "icarus_bnb_ccpi0_onbeam_syst_gundaminput.root"; // data = on-beam
  string outfile_signal_string = "icarus_bnb_ccpi0_signal_syst_gundaminput.root";

  // Test
  TChain ch("SelectedEvents");
  ch.Add((infile_string + selnutree_string).c_str());
  ch.Add((infile_string + selcostree_string).c_str());
  ch.Add((infile_string + seloffbeamtree_string).c_str());
  ROOT::RDataFrame rdf_mc(ch);
  ROOT::RDataFrame rdf_data(selonbeamtree_string, infile_string);
  ROOT::RDataFrame rdf_signal(signaltree_string, infile_string);

  // Recast "double" columns as "ints" and write to new TTree
  auto rdf_mc_conv = rdf_mc
    .Redefine("CutType", [](double d) { return static_cast<int>(d); }, {"CutType"})
    .Redefine("IsData", [](double d) { return static_cast<int>(d); }, {"IsData"})
    .Redefine("IsNu", [](double d) { return static_cast<int>(d); }, {"IsNu"})
    .Redefine("valid_syst", [](bool b) { return static_cast<int>(b); }, {"valid_syst"})
    .Redefine("category", [](double d) { return static_cast<int>(d); }, {"category"})
    .Redefine("category", "if (category < 0) return 8; else return category;") // 8 is cosmic category specified in spine_anaplot_tools
    .Filter("(IsNu == 1 && valid_syst == 1) || IsNu != 1") // remove duplicates that appear after ./run_systematics...
    .Snapshot("SelectedEvents", outfile_mc_string.c_str());
  
  auto rdf_data_conv = rdf_data
    .Redefine("CutType", [](double d) { return static_cast<int>(d); }, {"CutType"})
    .Redefine("IsData", [](double d) { return static_cast<int>(d); }, {"IsData"})
    .Redefine("IsNu", [](double d) { return static_cast<int>(d); }, {"IsNu"})
    .Redefine("valid_syst", [](bool b) { return static_cast<int>(b); }, {"valid_syst"})
    .Redefine("category", [](double d) { return static_cast<int>(d); }, {"category"})
    .Redefine("category", "if (category < 0) return 8; else return category;") // 8 is cosmic category specified in spine_anaplot_tools
    .Snapshot("SelectedEvents", outfile_data_string.c_str());

  auto rdf_signal_conv = rdf_signal
    .Redefine("CutType", [](double d) { return static_cast<int>(d); }, {"CutType"})
    .Redefine("IsData", [](double d) { return static_cast<int>(d); }, {"IsData"})
    .Redefine("IsNu", [](double d) { return static_cast<int>(d); }, {"IsNu"})
    .Redefine("valid_syst", [](bool b) { return static_cast<int>(b); }, {"valid_syst"})
    .Redefine("category", [](double d) { return static_cast<int>(d); }, {"category"})
    .Redefine("category", "if (category < 0) return 8; else return category;")
    .Filter("(IsNu == 1 && valid_syst == 1) || IsNu != 1")
    .Snapshot("SignalEvents", outfile_signal_string.c_str());

  // Write POT and Livetime of original samples to output file
  std::unique_ptr<TFile> infile( TFile::Open(infile_string.c_str(), "READ") );
  std::unique_ptr<TFile> outfile_mc( TFile::Open(outfile_mc_string.c_str(), "UPDATE") ); // "UPDATE" because we created these files above
  std::unique_ptr<TFile> outfile_data( TFile::Open(outfile_data_string.c_str(), "UPDATE") );

  TH1D *POT_mc = (TH1D*)infile->Get(potmc_string.c_str());
  TH1D *POT_mc_clone = (TH1D*)POT_mc->Clone("POT_mc");
  TH1D *POT_onbeam = (TH1D*)infile->Get(potonbeam_string.c_str());
  TH1D *POT_onbeam_clone = (TH1D*)POT_onbeam->Clone("POT_onbeam");
  TH1D *Livetime_mc = (TH1D*)infile->Get(livetimemc_string.c_str());
  TH1D *Livetime_mc_clone = (TH1D*)Livetime_mc->Clone("Livetime_mc");
  TH1D *Livetime_offbeam = (TH1D*)infile->Get(livetimeoffbeam_string.c_str());
  TH1D *Livetime_offbeam_clone = (TH1D*)Livetime_offbeam->Clone("Livetime_offbeam");
  TH1D *Livetime_onbeam = (TH1D*)infile->Get(livetimeonbeam_string.c_str());
  TH1D *Livetime_onbeam_clone = (TH1D*)Livetime_onbeam->Clone("Livetime_onbeam");

  outfile_mc->cd();
  POT_mc_clone->Write();
  Livetime_mc_clone->Write();
  Livetime_offbeam_clone->Write();

  outfile_data->cd();
  POT_onbeam_clone->Write();
  Livetime_onbeam_clone->Write();
  
}

vector<string> GetGENIEMultisigmaKnobNames()
{
  return 
    {
        "ZExpA1CCQE",
        "ZExpA2CCQE",
	"ZExpA3CCQE",
	"ZExpA4CCQE",
	"VecFFCCQEshape",
	"RPA_CCQE",
	 "CoulombCCQE",
	"NormCCMEC",
	"NormNCMEC",
	  //"DecayAngMEC",
	"MaNCEL",
	"EtaNCEL",
	"MaCCRES",
	"MvCCRES",
	"MaNCRES",
	"MvNCRES",
	"NonRESBGvpCC1pi",
	"NonRESBGvpCC2pi",
	"NonRESBGvpNC1pi",
	"NonRESBGvpNC2pi",
	"NonRESBGvnCC1pi",
	"NonRESBGvnCC2pi",
	"NonRESBGvnNC1pi" 
	  };
  
}
