#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

void treeReader(){
  TFile* f = TFile::Open("SLntuple.root");
  TIter nkey(f->GetListOfKeys());
  TKey* key = (TKey*)nkey();
  key->GetClassName();
  cout << "Directory key by name: " << key->GetName() << endl;;
  TDirectory* theDir = (TDirectory*)f->Get(key->GetName());
  TIter dkey(theDir->GetListOfKeys());
  TKey* ikey = (TKey*)dkey();
  ikey->GetClassName();
  cout << "Tree key by name: " << ikey->GetName() << endl;;
  TTree* tree = (TTree*)theDir->Get(ikey->GetName());
  tree->Scan();
}
