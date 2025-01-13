void sbruce2gundam()
{
  // Input/Output
  TFile* infile = new TFile("pi02024_output.root","READ");
  
  // Loop over TTree
  TTreeReader myReader("SelectedEvents_NSigmas", infile);
  TTreeReaderValue<Double_t> cut_type(myReader, "CutType");
  TTreeReaderValue<Double_t> is_signal(myReader, "IsSignal");
  TTreeReaderValue<Double_t> is_data(myReader, "IsData");
  ROOT::RVec<int> cuttype;
  ROOT::RVec<int> issignal;
  ROOT::RVec<int> isdata;
  while (myReader.Next())
    {

      // Convert double to int
      Int_t _cut_type = (Int_t) *cut_type;
      Int_t _is_signal = (Int_t) *is_signal;
      Int_t _is_data = (Int_t) *is_data;

      cuttype.push_back(_cut_type);
      issignal.push_back(_is_signal);
      isdata.push_back(_is_data);
    }

  // Replace (double) branch with new (int) branch
  ROOT::RDataFrame rdf("SelectedEvents_NSigmas", infile);
  rdf.Redefine("CutType", [cuttype](ULong64_t idx){return cuttype.at(idx);}, {"rdfentry_"})
     .Redefine("IsSignal", [issignal](ULong64_t idx){ return issignal.at(idx); }, {"rdfentry_"})
     .Redefine("IsData", [isdata](ULong64_t idx){ return isdata.at(idx); }, {"rdfentry_"})
     .Snapshot("SelectedEvents_NSigmas", "pi02024_output1.root");

  infile->Close();
}
