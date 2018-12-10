{
  //current rootlogon file. This one defines my top candidate packaging, a vector of pairs of TLorentzVectors and integers to be used as flags
  gSystem->Load("libFWCoreFWLite.so");
  gSystem->Load("libTreeViewer.so");
  //AutoLibraryLoader::enable();
  FWLiteEnabler::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  //  gROOT->ProcessLine(".include /cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/roofit/5.34.18-cms3/include/");
  gROOT->SetStyle ("Plain");
  gSystem->Load("libRooFit");
  //gInterpreter->GenerateDictionary("vector<Track>","Track.h;vector");
  gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
  gInterpreter->GenerateDictionary("vector<pair<vector<TLorentzVector>,vector<int>>>","TLorentzVector.h;vector;utility");
  using namespace RooFit ;
  cout << "loaded" << endl;
}
