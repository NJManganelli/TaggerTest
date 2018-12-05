#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"

//manditory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
//this include is useful to get the helper function to make the vector of constituents
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"

//this include is necessary to handle exceptions thrown by the top tagger code
#include "TopTagger/CfgParser/include/TTException.h"

#include "rootdict.h"

class ResTTEvaluator{
 public:
  ResTTEvaluator(std::string topTaggerName);
  ~ResTTEvaluator();

  //candidate setting and getting
  void addCand(std::vector<TLorentzVector>* cand, double discriminant);
  std::vector<  std::pair< std::vector <TLorentzVector>, double > > getAllCand();
  std::pair< std::vector <TLorentzVector>, double > getCand(int index);
  double getDiscr(int index);
  std::pair< std::vector <TLorentzVector>, double > getOrderedCand(int index);
  double getOrderedDiscr(int index);

  //Gen getting and setting //CHANGE: Make Gen strictly gen particles, Reco the Reco-matched gen particles
  //for pure Gen, do deltaR matching and permit boosted quark pairs to be matched, regardless of AK4 or AK8
  void addReco(std::vector<TLorentzVector>* gen);
  void addReco(std::vector<TLorentzVector>* gen, std::vector<int> flagVector); //using boost::(type)_bitset would be more space efficient ... up to factor of 8
  std::pair< std::vector<TLorentzVector>, std::vector<int> > getReco(int index);
  void addAllReco( std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > *allGenInput);

  void addGen(std::vector<TLorentzVector>* gen);
  void addGen(std::vector<TLorentzVector>* gen, std::vector<int> flagVector); //using boost::(type)_bitset would be more space efficient ... up to factor of 8
  std::pair< std::vector<TLorentzVector>, std::vector<int> > getGen(int index);
  void addAllGen( std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > *allGenInput);

  //print information, primarily for debugging at time of writing this comment
  void printCand(int index);
  void printGen(int index);
  void printMatrixR(int index); //the matching matrix showing which resolved top's gen-matched reco jets are the exact match for the candidate's constituents
  void printMatrixG(int index); //the matching matrix showing which resolved top's gen jets are the DeltaR match for the candidate's constituents
  void printDimensions(); //method to check internal dimensions, such as nReco, nGen, nCand, maxSize, match, etc.

  //Once all candidates have been added, this method evaluates everything
  void evaluateR();
  void evaluateG();
  
  
 private:
  const static uint _maxSize = 100;
  uint _nReco;
  uint _nGen;
  uint _nCand;
  uint _matchR[_maxSize][_maxSize][3]; //support max 100 candidates and gen tops... shouldn't ever be a problem at 13TeV CoM Energy 
  uint _matchG[_maxSize][_maxSize][3]; //support max 100 candidates and gen tops... shouldn't ever be a problem at 13TeV CoM Energy 
  //encode in this matrix the info for whether it was even POSSIBLE to match that jet... 
  bool _haveEvaluatedR;
  bool _haveEvaluatedG;
  std::string _topTaggerName;
  uint _orderIndex[_maxSize];
  std::vector< std::pair< std::vector<TLorentzVector>, double > > _cand;
  std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > _reco; //include flags in here...
  std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > _gen; //include flags in here...
  //std::vector<std::vector<uint>> _flags; //Have to decide how to order these things in general way...
  bool _haveFlags;
  std::vector<std::string> _matchClassification;
};
ResTTEvaluator::ResTTEvaluator(std::string topTaggerName){
  _nReco = 0;
  _nGen = 0;
  _nCand = 0;
  _haveFlags = false;
  _topTaggerName = topTaggerName;
  std::cout << "ResTTEvaluator being created" << std::endl;
}
ResTTEvaluator::~ResTTEvaluator(){
  std::cout << "ResTTEvaluator being destroyed" << std::endl;
}
void ResTTEvaluator::addCand(std::vector<TLorentzVector>* cand, double discriminant){
  std::pair< std::vector<TLorentzVector>, double> temp;
  temp.first = *cand;
  temp.second = discriminant;
  _cand.push_back(temp);
  _nCand++;
}
std::vector<std::pair<std::vector<TLorentzVector>, double>> ResTTEvaluator::getAllCand(){
  return _cand;
}
std::pair< std::vector<TLorentzVector>, double > ResTTEvaluator::getCand(int index){
  return _cand[index];
}
double ResTTEvaluator::getDiscr(int index){
  return _cand[index].second;
}
std::pair< std::vector<TLorentzVector>, double > ResTTEvaluator::getOrderedCand(int index){
  std::cout << "Method not implemented, just returning unordered candidate..." << std::endl;
  return _cand[index];
}
double ResTTEvaluator::getOrderedDiscr(int index){
  std::cout << "Method not implemented, just returning unordered candidate's discriminant..." << std::endl;
  return _cand[index].second;
}
void ResTTEvaluator::addReco(std::vector<TLorentzVector>* genReco){
  if(_haveFlags){
    _haveFlags = false;
    std::cout << "Added genReco without flags, but previous candidate had flags! Will not perform any calculations depending on presence of flags..." << std::endl;
  }
  std::vector<int> tempFlags;
  tempFlags.push_back(-1);
  std::pair<std::vector<TLorentzVector>, std::vector<int> > tempReco;
  tempReco.first = *genReco;
  tempReco.second = tempFlags;
  _reco.push_back(tempReco);
  _nReco++;
}
void ResTTEvaluator::addReco(std::vector<TLorentzVector>* genReco, std::vector<int> flags){
  if(_reco.size() > 0 && !_haveFlags){
    _haveFlags = false; //explicit for readibility
    std::cout << "Added genReco with flags, but previous candidate didn't have flags! Will not perform any calculations depending on presence of flags..." << std::endl;
  }
  else
    _haveFlags = true;
  std::pair<std::vector<TLorentzVector>, std::vector<int> > tempReco;
  tempReco.first = *genReco;
  tempReco.second = flags;
  _reco.push_back(tempReco);
  _nReco++;
}
void ResTTEvaluator::addAllReco(std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > *allRecoInput){
  if(_reco.size() == 0){
    _reco = *allRecoInput;
    _nReco = _reco.size();
  }
  else{
    std::cout << "Attempted to add all genRecos when some have already been added. I can't let you do that, Dave" << std::endl;
    _haveFlags = true;
  }
}
std::pair< std::vector<TLorentzVector>, std::vector<int> > ResTTEvaluator::getReco(int index){
  return _reco[index];
}
void ResTTEvaluator::addGen(std::vector<TLorentzVector>* gen){
  if(_haveFlags){
    _haveFlags = false;
    std::cout << "Added gen without flags, but previous candidate had flags! Will not perform any calculations depending on presence of flags..." << std::endl;
  }
  else
    _haveFlags = true;
  std::vector<int> tempFlags;
  tempFlags.push_back(-1);
  std::pair<std::vector<TLorentzVector>, std::vector<int> > tempGen;
  tempGen.first = *gen;
  tempGen.second = tempFlags;
  _gen.push_back(tempGen);
  _nGen++;
}
void ResTTEvaluator::addGen(std::vector<TLorentzVector>* gen, std::vector<int> flags){
  if(_gen.size() > 0 && !_haveFlags){
    _haveFlags = false; //explicit for readibility
    std::cout << "Added gen with flags, but previous candidate didn't have flags! Will not perform any calculations depending on presence of flags..." << std::endl;
  }
  else
    _haveFlags = true;
  std::pair<std::vector<TLorentzVector>, std::vector<int> > tempGen;
  tempGen.first = *gen;
  tempGen.second = flags;
  _gen.push_back(tempGen);
  _nGen++;
}
void ResTTEvaluator::addAllGen(std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > *allGenInput){
  if(_gen.size() == 0){
    _gen = *allGenInput;
    _nGen = _gen.size();
  }
  else{
    std::cout << "Attempted to add all gens when some have already been added. I can't let you do that, Dave" << std::endl;
    _haveFlags = true;
  }
}
std::pair< std::vector<TLorentzVector>, std::vector<int> > ResTTEvaluator::getGen(int index){
  return _gen[index];
}
void ResTTEvaluator::printCand(int index){
  std::cout << "This isn't implemented yet... " << std::endl;
}
void ResTTEvaluator::printGen(int index){
  std::cout << "This isn't implemented yet... " << std::endl;
}
void ResTTEvaluator::printMatrixR(int index){
  if(_nReco > 0){
    std::cout << "\n\t Match Matrix for Top Candidate " << index << "\n\t\t\tb jet\tq1 jet\tq2 jet";
    for(int i = 0; i < _nReco; i++)
      std::cout << "\n\t Gen Reco " << i+1 << "\t" << _matchR[index][i][0] << "\t" << _matchR[index][i][1] << "\t" << _matchR[index][i][2];
    std::cout << std::endl;
  }
}
void ResTTEvaluator::printMatrixG(int index){
  if(_nGen > 0){
    std::cout << "\n\t Match Matrix for Top Candidate " << index << "\n\t\t\tb jet\tq1 jet\tq2 jet";
    for(int i = 0; i < _nGen; i++)
      std::cout << "\n\t Gen Part " << i+1 << "\t" << _matchG[index][i][0] << "\t" << _matchG[index][i][1] << "\t" << _matchG[index][i][2];
    std::cout << std::endl;
  }
}
void ResTTEvaluator::printDimensions(){
  std::cout << "nGen: " << _nGen << "\tnReco: " << _nReco << "\tnCand: " << _nCand << "\tmaxSize: " << _maxSize << std::endl;
}
void ResTTEvaluator::evaluateR(){
  for(int m = 0; m < _nCand; m++) //evaluate for each candidate
    for(int n = 0; n < _cand[m].first.size(); n++) //each candidate's constituent
      //for(int ii = 0; ii < 3; ii++) //each candidate's constituent
      for(int o = 0; o < _nReco; o++) //for each candidate, check each Reconstructed top quark's constituents
	for(int p = 0; p < _reco[o].first.size(); p++) //check all the jets in the collection... should be 3 at all times... need to protect...
	//for(int jj = 0; jj < 3; jj++) //check all the jets in the collection... should be 3 at all times... need to protect...
	  {
	    std::cout << "\n\ti|ii|j|jj: " << m << n << o << p << "\t" << _cand[m].first[n].Pt() << "\t" << _reco[o].first[p].Pt();
	    _matchR[m][o][p] = (_reco[o].first[p] == _cand[m].first[n]); //for each candidate, store match truth in the matrix
	    std::cout << "\t" << _matchR[m][o][p];

	  }
  _haveEvaluatedR = true;
  std::cout << std::endl;
}
void ResTTEvaluator::evaluateG(){
  std::cout << "If I were a real little evaluator (Reco), I would have done some evaluation (and stuff!). Since I'm not (yet), I'll just let you know this worked!" << std::endl;
}

int main()
{
  //bool for silencing original top quark properties (except event #)
  bool silentRunning = true;
  bool debug1 = false;
  bool debug2 = false;
  //silentRunning = true;
  bool runStdExample = false;
  TFile *tf, *tf2, *of;
  TDirectory *td;
  TTree *tree;
  //std::string postfix = "tttt";

  //HOT discriminant histos
  TH1F *h_typeIII_hot = new TH1F ("h_typeIII_hot_", "Type III (correct) Top Quarks; Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0); 
  TH1F *h_typeIIb_hot = new TH1F ("h_typeIIb_hot_", "Type II (b swapped) Top Quarks; Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0); 
  TH1F *h_typeIIw_hot = new TH1F ("h_typeIIw_hot_", "Type II (q1 or q2 swapped) Top Quarks; Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0); 
  TH1F *h_typeIIo_hot = new TH1F ("h_typeIIo_hot_", "Type II (other) Top Quarks; Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0); 
  TH1F *h_typeI_hot = new TH1F ("h_typeI_hot_", "Type I (2+ top-daughters matched, 1 per reco top); Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_type0_hot = new TH1F ("h_type0_hot_", "Type 0 (all other tagger candidates); Discriminant ; Number of Tagger Candidates", 20, 0.0, 1.0);

  TH1F *h_eventN1_III_hot = new TH1F ("h_eventN1_tIII_hot", "Highest Disc Cand (Type III) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN1_II_hot = new TH1F ("h_eventN1_tII_hot", "Highest Disc Cand (Type II) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN1_I_hot = new TH1F ("h_eventN1_tI_hot", "Highest Disc Cand (Type I) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN1_0_hot = new TH1F ("h_eventN1_t0_hot", "Highest Disc Cand (Type 0) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);

  TH1F *h_eventN2_III_hot = new TH1F ("h_eventN2_tIII_hot", "2nd Highest Disc Cand (Type III) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN2_II_hot = new TH1F ("h_eventN2_tII_hot", "2nd Highest Disc Cand (Type II) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN2_I_hot = new TH1F ("h_eventN2_tI_hot", "2nd Highest Disc Cand (Type I) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN2_0_hot = new TH1F ("h_eventN2_t0_hot", "2nd Highest Disc Cand (Type 0) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);

  TH1F *h_eventN3_III_hot = new TH1F ("h_eventN3_tIII_hot", "3rd Highest Disc Cand (Type III) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN3_II_hot = new TH1F ("h_eventN3_tII_hot", "3rd Highest Disc Cand (Type II) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN3_I_hot = new TH1F ("h_eventN3_tI_hot", "3rd Highest Disc Cand (Type I) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN3_0_hot = new TH1F ("h_eventN3_t0_hot", "3rd Highest Disc Cand (Type 0) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);

  //BDT Candidates
  TH1F *h_typeIII_bdt = new TH1F ("h_typeIII_bdt_", "Type III (correct) Top Quarks; Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0); 
  TH1F *h_typeIIb_bdt = new TH1F ("h_typeIIb_bdt_", "Type II (b swapped) Top Quarks; Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0); 
  TH1F *h_typeIIw_bdt = new TH1F ("h_typeIIw_bdt_", "Type II (q1 or q2 swapped) Top Quarks; Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0); 
  TH1F *h_typeIIo_bdt = new TH1F ("h_typeIIo_bdt_", "Type II (other) Top Quarks; Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0); 
  TH1F *h_typeI_bdt = new TH1F ("h_typeI_bdt_", "Type I (2+ top-daughters matched, 1 per reco top); Discriminant; Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_type0_bdt = new TH1F ("h_type0_bdt_", "Type 0 (all other tagger candidates); Discriminant ; Number of Tagger Candidates", 20, 0.0, 1.0);

  TH1F *h_eventN1_III_bdt = new TH1F ("h_eventN1_tIII_bdt", "Highest Disc Cand (Type III) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN1_II_bdt = new TH1F ("h_eventN1_tII_bdt", "Highest Disc Cand (Type II) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN1_I_bdt = new TH1F ("h_eventN1_tI_bdt", "Highest Disc Cand (Type I) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN1_0_bdt = new TH1F ("h_eventN1_t0_bdt", "Highest Disc Cand (Type 0) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);

  TH1F *h_eventN2_III_bdt = new TH1F ("h_eventN2_tIII_bdt", "2nd Highest Disc Cand (Type III) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN2_II_bdt = new TH1F ("h_eventN2_tII_bdt", "2nd Highest Disc Cand (Type II) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN2_I_bdt = new TH1F ("h_eventN2_tI_bdt", "2nd Highest Disc Cand (Type I) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN2_0_bdt = new TH1F ("h_eventN2_t0_bdt", "2nd Highest Disc Cand (Type 0) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);

  TH1F *h_eventN3_III_bdt = new TH1F ("h_eventN3_tIII_bdt", "3rd Highest Disc Cand (Type III) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN3_II_bdt = new TH1F ("h_eventN3_tII_bdt", "3rd Highest Disc Cand (Type II) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN3_I_bdt = new TH1F ("h_eventN3_tI_bdt", "3rd Highest Disc Cand (Type I) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);
  TH1F *h_eventN3_0_bdt = new TH1F ("h_eventN3_t0_bdt", "3rd Highest Disc Cand (Type 0) ;Discriminant; Number of Number of Tagger Candidates", 20, 0.0, 1.0);  

  //TH1F *TypeWDiscr = new TH1F ("h_fullWdiscr", "Fully Reconstructible W bosons; Discriminant; Number of matches", 20, 0.0, 1.0);

  TH1I *RecoTypes = new TH1I ("h_recoTypes", "Reconstructible Top Candidates; I: IIbq: IIqq: III:;Number Gen Particles", 3, 0, 3);
  std::cout << "Testing the file for Top Quarks" << std::endl;
  
  if(runStdExample == true){
    //Open input ntuple file 
    tf = TFile::Open("exampleInputs.root");

    //Get tree from file
    tree = (TTree*)tf->Get("slimmedTuple");
  }
  else if(runStdExample == false){
    //Open input ntuple file

    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/SLntuple.root");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/SLntupleMasked.root");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/200kFourTop.root");
    tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/200kTwoTop.root");

    //Get TDirectory next
    td = (TDirectory*)tf2->Get("tree");

    //Get tree from file
    tree = (TTree*)td->Get("nTuple");
  }
  std::cout << "Loaded the file" << std::endl;


    //Variables to hold inputs
    //AK4 jet variables
    //Each entry in these vectors refers to information for 1 AK4 jet
    std::vector<TLorentzVector>** AK4JetLV = new std::vector<TLorentzVector>*();
    std::vector<float>** AK4JetBtag = new std::vector<float>*();
    std::vector<float>** AK4qgMult = new std::vector<float>*();
    std::vector<float>** AK4qgPtD = new std::vector<float>*();
    std::vector<float>** AK4qgAxis1 = new std::vector<float>*();
    std::vector<float>** AK4qgAxis2 = new std::vector<float>*();
    std::vector<float>** AK4recoJetschargedHadronEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4recoJetschargedEmEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4recoJetsneutralEmEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4ElectronEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4PhotonEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4recoJetsneutralEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4recoJetsHFHadronEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4recoJetsmuonEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4recoJetsHFEMEnergyFraction = new std::vector<float>*();
    std::vector<float>** AK4NeutralHadronMultiplicity = new std::vector<float>*();
    std::vector<float>** AK4ChargedHadronMultiplicity = new std::vector<float>*();
    std::vector<float>** AK4ElectronMultiplicity = new std::vector<float>*();
    std::vector<float>** AK4MuonMultiplicity = new std::vector<float>*();
    std::vector<float>** AK4PhotonMultiplicity = new std::vector<float>*();
    std::vector<float>** AK4DeepCSVbb = new std::vector<float>*();
    std::vector<float>** AK4DeepCSVb = new std::vector<float>*();
    std::vector<float>** AK4DeepCSVc = new std::vector<float>*();
    std::vector<float>** AK4DeepCSVcc = new std::vector<float>*();
    std::vector<float>** AK4DeepCSVl = new std::vector<float>*();
    std::cout << "Created all AK4 jet variables" << std::endl;
    
    //AK8 jet varaibles
    //The elements of each vector refer to one AK8 jet
    std::vector<TLorentzVector>** AK8JetLV = new std::vector<TLorentzVector>*();
    std::vector<float>** AK8JetSoftdropMass = new std::vector<float>*();
    std::vector<float>** AK8JetDeepAK8Top = new std::vector<float>*();
    std::vector<float>** AK8JetDeepAK8W = new std::vector<float>*();
    std::vector<std::vector<TLorentzVector>>** AK8SubjetLV = new std::vector<std::vector<TLorentzVector>>*();
    std::cout << "Created all AK8 jet variables, preparing to Deactivate then reactivate branches..." << std::endl;

    //Create custom variables for reading in top gen-matched jets and status flags
    std::vector<TLorentzVector>** hadTop1Constit = new std::vector<TLorentzVector>*();
    std::vector<TLorentzVector>** hadTop2Constit = new std::vector<TLorentzVector>*();
    std::vector<TLorentzVector>** hadTop3Constit = new std::vector<TLorentzVector>*();
    //std::vector<TLorentzVector>** lepTopConstit = new std::vector<TLorentzVector>*(); //for future improvement
    std::vector<uint>** FlagTop = new std::vector<uint>*();
    std::vector<uint>** FlagBottom = new std::vector<uint>*();
    std::vector<uint>** FlagQ1 = new std::vector<uint>*();
    std::vector<uint>** FlagQ2 = new std::vector<uint>*();

    //Deactivate all branches, then activate the branches of interest
    tree->SetBranchStatus("*", 0);
    std::cout << "Deactivated all branches" << std::endl;

    //standard exampleInputs.root settings
    if(runStdExample == true){
      //Activate branches of interest
      //AK4 jet lorentz vectors
      tree->SetBranchStatus( "ak4jetsLVec", 1);
      tree->SetBranchAddress("ak4jetsLVec", AK4JetLV);
    
      //AK4 jet b-tag values (0 not a b, 1 is a b)
      tree->SetBranchStatus( "ak4recoJetsBtag", 1);
      tree->SetBranchAddress("ak4recoJetsBtag", AK4JetBtag);

      //AK4 qg jet multiplicity
      tree->SetBranchStatus( "ak4qgMult", 1);
      tree->SetBranchAddress("ak4qgMult", AK4qgMult);

      //AK4 qg PtD
      tree->SetBranchStatus( "ak4qgPtD", 1);
      tree->SetBranchAddress("ak4qgPtD", AK4qgPtD);

      //AK4 qg jet semimajor axis
      tree->SetBranchStatus( "ak4qgAxis1", 1);
      tree->SetBranchAddress("ak4qgAxis1", AK4qgAxis1);

      //AK4 qg jet semiminor axis
      tree->SetBranchStatus( "ak4qgAxis2", 1);
      tree->SetBranchAddress("ak4qgAxis2", AK4qgAxis2);

      //AK4 jet charged hadronic energy fraction 
      tree->SetBranchStatus( "ak4recoJetschargedHadronEnergyFraction", 1);
      tree->SetBranchAddress("ak4recoJetschargedHadronEnergyFraction", AK4recoJetschargedHadronEnergyFraction);

      //AK4 jet charged electromagnetic energy fraction 
      tree->SetBranchStatus( "ak4recoJetschargedEmEnergyFraction", 1);
      tree->SetBranchAddress("ak4recoJetschargedEmEnergyFraction", AK4recoJetschargedEmEnergyFraction);

      //AK4 jet neutral hadronic energy fraction 
      tree->SetBranchStatus( "ak4recoJetsneutralEmEnergyFraction", 1);
      tree->SetBranchAddress("ak4recoJetsneutralEmEnergyFraction", AK4recoJetsneutralEmEnergyFraction);

      //AK4 electron energy fraction
      tree->SetBranchStatus( "ak4ElectronEnergyFraction", 1);
      tree->SetBranchAddress("ak4ElectronEnergyFraction", AK4ElectronEnergyFraction);

      //AK4 photon energy fraction 
      tree->SetBranchStatus( "ak4PhotonEnergyFraction", 1);
      tree->SetBranchAddress("ak4PhotonEnergyFraction", AK4PhotonEnergyFraction);

      //AK4 neutral energy fraction
      tree->SetBranchStatus( "ak4recoJetsneutralEnergyFraction", 1);
      tree->SetBranchAddress("ak4recoJetsneutralEnergyFraction", AK4recoJetsneutralEnergyFraction);

      //AK4 HF hadronic energy fraction 
      tree->SetBranchStatus( "ak4recoJetsHFHadronEnergyFraction", 1);
      tree->SetBranchAddress("ak4recoJetsHFHadronEnergyFraction", AK4recoJetsHFHadronEnergyFraction);

      //AK4 muon energy fraction 
      tree->SetBranchStatus( "ak4recoJetsmuonEnergyFraction", 1);
      tree->SetBranchAddress("ak4recoJetsmuonEnergyFraction", AK4recoJetsmuonEnergyFraction);

      //AK4 HF electron energy fraction 
      tree->SetBranchStatus( "ak4recoJetsHFEMEnergyFraction", 1);
      tree->SetBranchAddress("ak4recoJetsHFEMEnergyFraction", AK4recoJetsHFEMEnergyFraction);

      //AK4 neutral hadronic particle multiplicity
      tree->SetBranchStatus( "ak4NeutralHadronMultiplicity", 1);
      tree->SetBranchAddress("ak4NeutralHadronMultiplicity", AK4NeutralHadronMultiplicity);

      //AK4 charged hadronic particle multiplicity
      tree->SetBranchStatus( "ak4ChargedHadronMultiplicity", 1);
      tree->SetBranchAddress("ak4ChargedHadronMultiplicity", AK4ChargedHadronMultiplicity);

      //AK4 electron multiplicity
      tree->SetBranchStatus( "ak4ElectronMultiplicity", 1);
      tree->SetBranchAddress("ak4ElectronMultiplicity", AK4ElectronMultiplicity);

      //AK4 muon multiplicity
      tree->SetBranchStatus( "ak4MuonMultiplicity", 1);
      tree->SetBranchAddress("ak4MuonMultiplicity", AK4MuonMultiplicity);

      //AK4 photon multiplicity
      tree->SetBranchStatus( "ak4PhotonMultiplicity", 1);
      tree->SetBranchAddress("ak4PhotonMultiplicity", AK4PhotonMultiplicity);

      //AK4 deepCSV bb discriminator
      tree->SetBranchStatus( "ak4DeepCSVbb", 1);
      tree->SetBranchAddress("ak4DeepCSVbb", AK4DeepCSVbb);

      //AK4 deepCSV b discriminator
      tree->SetBranchStatus( "ak4DeepCSVb", 1);
      tree->SetBranchAddress("ak4DeepCSVb", AK4DeepCSVb);

      //AK4 deepCSV charm discriminator
      tree->SetBranchStatus( "ak4DeepCSVc", 1);
      tree->SetBranchAddress("ak4DeepCSVc", AK4DeepCSVc);

      //AK4 deepCSV double charm discriminator
      tree->SetBranchStatus( "ak4DeepCSVcc", 1);
      tree->SetBranchAddress("ak4DeepCSVcc", AK4DeepCSVcc);

      //AK4 deepCSV light discriminator
      tree->SetBranchStatus( "ak4DeepCSVl", 1);
      tree->SetBranchAddress("ak4DeepCSVl", AK4DeepCSVl);

      //HIDE THESE and check    
      // //AK8 jet lorentz vectors
      // tree->SetBranchStatus( "ak8JetsLVec", 0);
      // tree->SetBranchAddress("ak8JetsLVec", AK8JetLV);
    
      // //AK8 subjet lorentz vectors (soft drop algo produces 2 subjets for each AK8 jet)
      // tree->SetBranchStatus( "ak8SubJetsLVec", 0);
      // tree->SetBranchAddress("ak8SubJetsLVec", AK8SubjetLV);

      // //AK8 jet deepAK8 top discriminator variable
      // tree->SetBranchStatus( "ak8DeepAK8Top", 0);
      // tree->SetBranchAddress("ak8DeepAK8Top", AK8JetDeepAK8Top);
    
      // //AK8 jet deepAK8 W discriminator variable
      // tree->SetBranchStatus( "ak8DeepAK8W", 0);
      // tree->SetBranchAddress("ak8DeepAK8W", AK8JetDeepAK8W);
    
      // //AK8 jet softdrop mass
      // tree->SetBranchStatus( "ak8softDropMass", 0);
      // tree->SetBranchAddress("ak8softDropMass", AK8JetSoftdropMass);
    }
    else if(runStdExample == false){
      //SLntuple.root branches
      //Activate branches of interest
      //AK4 jet lorentz vectors
      tree->SetBranchStatus( "JetLVec", 1);
      tree->SetBranchAddress("JetLVec", AK4JetLV);
    
      //AK4 jet b-tag values (0 not a b, 1 is a b)
      tree->SetBranchStatus( "btag", 1);
      tree->SetBranchAddress("btag", AK4JetBtag);

      //AK4 qg jet multiplicity
      tree->SetBranchStatus( "qgMult", 1);
      tree->SetBranchAddress("qgMult", AK4qgMult);

      //AK4 qg PtD
      tree->SetBranchStatus( "qgPtD", 1);
      tree->SetBranchAddress("qgPtD", AK4qgPtD);

      //AK4 qg jet semimajor axis
      tree->SetBranchStatus( "qgAxis1", 1);
      tree->SetBranchAddress("qgAxis1", AK4qgAxis1);

      //AK4 qg jet semiminor axis
      tree->SetBranchStatus( "qgAxis2", 1);
      tree->SetBranchAddress("qgAxis2", AK4qgAxis2);

      //AK4 jet charged hadronic energy fraction 
      tree->SetBranchStatus( "chargedHadronEnergyFraction", 1);
      tree->SetBranchAddress("chargedHadronEnergyFraction", AK4recoJetschargedHadronEnergyFraction);

      //AK4 jet charged electromagnetic energy fraction 
      tree->SetBranchStatus( "chargedEmEnergyFraction", 1);
      tree->SetBranchAddress("chargedEmEnergyFraction", AK4recoJetschargedEmEnergyFraction);

      //AK4 jet neutral hadronic energy fraction 
      tree->SetBranchStatus( "neutralEmEnergyFraction", 1);
      tree->SetBranchAddress("neutralEmEnergyFraction", AK4recoJetsneutralEmEnergyFraction);

      //AK4 electron energy fraction
      tree->SetBranchStatus( "electronEnergyFraction", 1);
      tree->SetBranchAddress("electronEnergyFraction", AK4ElectronEnergyFraction);

      //AK4 photon energy fraction 
      tree->SetBranchStatus( "photonEnergyFraction", 1);
      tree->SetBranchAddress("photonEnergyFraction", AK4PhotonEnergyFraction);

      //AK4 neutral energy fraction #WHAT IS THIS
      tree->SetBranchStatus( "neutralHadronEnergyFraction", 1);
      tree->SetBranchAddress("neutralHadronEnergyFraction", AK4recoJetsneutralEnergyFraction);

      //AK4 HF hadronic energy fraction 
      tree->SetBranchStatus( "recoJetsHFHadronEnergyFraction", 1);
      tree->SetBranchAddress("recoJetsHFHadronEnergyFraction", AK4recoJetsHFHadronEnergyFraction);

      //AK4 muon energy fraction 
      tree->SetBranchStatus( "muonEnergyFraction", 1);
      tree->SetBranchAddress("muonEnergyFraction", AK4recoJetsmuonEnergyFraction);

      //AK4 HF electron energy fraction 
      tree->SetBranchStatus( "recoJetsHFEMEnergyFraction", 1);
      tree->SetBranchAddress("recoJetsHFEMEnergyFraction", AK4recoJetsHFEMEnergyFraction);

      //AK4 neutral hadronic particle multiplicity
      tree->SetBranchStatus( "neutralHadronMultiplicity", 1);
      tree->SetBranchAddress("neutralHadronMultiplicity", AK4NeutralHadronMultiplicity);

      //AK4 charged hadronic particle multiplicity
      tree->SetBranchStatus( "chargedHadronMultiplicity", 1);
      tree->SetBranchAddress("chargedHadronMultiplicity", AK4ChargedHadronMultiplicity);

      //AK4 electron multiplicity
      tree->SetBranchStatus( "electronMultiplicity", 1);
      tree->SetBranchAddress("electronMultiplicity", AK4ElectronMultiplicity);

      //AK4 muon multiplicity
      tree->SetBranchStatus( "muonMultiplicity", 1);
      tree->SetBranchAddress("muonMultiplicity", AK4MuonMultiplicity);

      //AK4 photon multiplicity
      tree->SetBranchStatus( "photonMultiplicity", 1);
      tree->SetBranchAddress("photonMultiplicity", AK4PhotonMultiplicity);

      //AK4 deepCSV bb discriminator
      tree->SetBranchStatus( "deepCSVbb", 1);
      tree->SetBranchAddress("deepCSVbb", AK4DeepCSVbb);

      //AK4 deepCSV b discriminator
      tree->SetBranchStatus( "deepCSVb", 1);
      tree->SetBranchAddress("deepCSVb", AK4DeepCSVb);

      //AK4 deepCSV charm discriminator
      tree->SetBranchStatus( "deepCSVc", 1);
      tree->SetBranchAddress("deepCSVc", AK4DeepCSVc);

      //AK4 deepCSV double charm discriminator
      tree->SetBranchStatus( "deepCSVcc", 1);
      tree->SetBranchAddress("deepCSVcc", AK4DeepCSVcc);

      //AK4 deepCSV light discriminator
      tree->SetBranchStatus( "deepCSVl", 1);
      tree->SetBranchAddress("deepCSVl", AK4DeepCSVl);

      //Hadronic Top candidate 1 vector of TLV's
      tree->SetBranchStatus( "hadTop1Constit", 1);
      tree->SetBranchAddress("hadTop1Constit", hadTop1Constit);

      //Hadronic Top candidate 2 vector of TLV's
      tree->SetBranchStatus( "hadTop2Constit", 1);
      tree->SetBranchAddress("hadTop2Constit", hadTop2Constit);

      //Hadronic Top candidate 3 vector of TLV's
      tree->SetBranchStatus( "hadTop3Constit", 1);
      tree->SetBranchAddress("hadTop3Constit", hadTop3Constit);

      //Top flags for overall object. >=1000 indicates leptonic top
      tree->SetBranchStatus( "FlagTop", 1);
      tree->SetBranchAddress("FlagTop", FlagTop);

      //Bottom flags, Reconstructible if >= ???
      tree->SetBranchStatus( "FlagBottom", 1);
      tree->SetBranchAddress("FlagBottom", FlagBottom);

      //Q1 flags, Reconstructible if >= ???
      tree->SetBranchStatus( "FlagQ1", 1);
      tree->SetBranchAddress("FlagQ1", FlagQ1);

      //Q2 flags, Reconstructible if >= ???
      tree->SetBranchStatus( "FlagQ2", 1);
      tree->SetBranchAddress("FlagQ2", FlagQ2);

      std::cout << "tree branches attached" << std::endl;
    }

    //Create top tagger object
    TopTagger tt;

    //try-catch on TTException which are thrown by the top tagger
    try
    {
        //Set top tagger cfg file
        tt.setCfgFile("TopTagger.cfg");

        //Loop over events
        int Nevt = 0;
        while(tree->GetEntry(Nevt))
        {
            //increment event number
            ++Nevt;

            //Print event number 
            printf("Event #: %i\n", Nevt);

            //Use helper function to create input list 
            //Create AK4 inputs object
            ttUtility::ConstAK4Inputs<float> AK4Inputs(**AK4JetLV, **AK4JetBtag);
            AK4Inputs.addSupplamentalVector("qgPtD",                                **AK4qgPtD);
            AK4Inputs.addSupplamentalVector("qgAxis1",                              **AK4qgAxis1);
            AK4Inputs.addSupplamentalVector("qgAxis2",                              **AK4qgAxis2);
            AK4Inputs.addSupplamentalVector("qgMult",                               **AK4qgMult);
            AK4Inputs.addSupplamentalVector("recoJetschargedHadronEnergyFraction",  **AK4recoJetschargedHadronEnergyFraction);
            AK4Inputs.addSupplamentalVector("recoJetschargedEmEnergyFraction",      **AK4recoJetschargedEmEnergyFraction);
            AK4Inputs.addSupplamentalVector("recoJetsneutralEmEnergyFraction",      **AK4recoJetsneutralEmEnergyFraction);
            AK4Inputs.addSupplamentalVector("recoJetsmuonEnergyFraction",           **AK4recoJetsmuonEnergyFraction);
            AK4Inputs.addSupplamentalVector("recoJetsHFHadronEnergyFraction",       **AK4recoJetsHFHadronEnergyFraction);
            AK4Inputs.addSupplamentalVector("recoJetsHFEMEnergyFraction",           **AK4recoJetsHFEMEnergyFraction);
            AK4Inputs.addSupplamentalVector("recoJetsneutralEnergyFraction",        **AK4recoJetsneutralEnergyFraction);
            AK4Inputs.addSupplamentalVector("PhotonEnergyFraction",                 **AK4PhotonEnergyFraction);
            AK4Inputs.addSupplamentalVector("ElectronEnergyFraction",               **AK4ElectronEnergyFraction);
            AK4Inputs.addSupplamentalVector("ChargedHadronMultiplicity",            **AK4ChargedHadronMultiplicity);
            AK4Inputs.addSupplamentalVector("NeutralHadronMultiplicity",            **AK4NeutralHadronMultiplicity);
            AK4Inputs.addSupplamentalVector("PhotonMultiplicity",                   **AK4PhotonMultiplicity);
            AK4Inputs.addSupplamentalVector("ElectronMultiplicity",                 **AK4ElectronMultiplicity);
            AK4Inputs.addSupplamentalVector("MuonMultiplicity",                     **AK4MuonMultiplicity);
            AK4Inputs.addSupplamentalVector("DeepCSVb",                             **AK4DeepCSVb);
            AK4Inputs.addSupplamentalVector("DeepCSVc",                             **AK4DeepCSVc);
            AK4Inputs.addSupplamentalVector("DeepCSVl",                             **AK4DeepCSVl);
            AK4Inputs.addSupplamentalVector("DeepCSVbb",                            **AK4DeepCSVbb);
            AK4Inputs.addSupplamentalVector("DeepCSVcc",                            **AK4DeepCSVcc);

            
            //Create AK8 inputs object
            // ttUtility::ConstAK8Inputs<float> AK8Inputs(
            //     **AK8JetLV,
            //     **AK8JetDeepAK8Top,
            //     **AK8JetDeepAK8W,
            //     **AK8JetSoftdropMass,
            //     **AK8SubjetLV
            //     );
            
            //Create jets constituents list combining AK4 and AK8 jets, these are used to construct top candiates
            //The vector of input constituents can also be constructed "by hand"
            // std::vector<Constituent> constituents = ttUtility::packageConstituents(AK4Inputs, AK8Inputs);
            std::vector<Constituent> constituents = ttUtility::packageConstituents(AK4Inputs);

            //run the top tagger
            tt.runTagger(std::move(constituents));

            //retrieve the top tagger results object
            const TopTaggerResults& ttr = tt.getResults();

            //get reconstructed top
            const std::vector<TopObject*>& tops = ttr.getTops();

            //print the number of tops found in the event 
            if(!silentRunning) printf("\tN tops: %ld\n", tops.size());
		
	    //count reconstructible tops via the flags...
	    uint nRecoTops = 0;
	    uint nRecoWs = 0;
	    bool Top1Reco = false;
	    bool Top2Reco = false;
	    bool Top3Reco = false;
	    if(debug1) printf("Booleans for Top reconstruction (initialization) %2d %2d %2d", Top1Reco, Top2Reco, Top3Reco);
	    uint flagcounter = 0;
	    if(debug1) printf("\n\tDebugging reconstructible top counting and flags: ");
	    ResTTEvaluator HOTEval("HOT");
	    HOTEval.printCand(1);
	    HOTEval.printGen(2);
	    std::vector<int> Top1flags;
	    Top1flags.push_back(1);
	    Top1flags.push_back(2);
	    HOTEval.addReco(*hadTop2Constit, Top1flags);
	    HOTEval.addReco(*hadTop1Constit, Top1flags);
	    HOTEval.addReco(*hadTop2Constit, Top1flags);
	    HOTEval.addReco(*hadTop3Constit, Top1flags);
	    HOTEval.addCand(*hadTop1Constit, 0.873);
	    HOTEval.addCand(*hadTop2Constit, 0.997);
	    HOTEval.addCand(*hadTop3Constit, 0.469);
	    HOTEval.addCand(*hadTop2Constit, 0.371);
	    HOTEval.addCand(*hadTop1Constit, 0.214);
	    std::cout << "\n=======================================" << std::endl;
	    if((*hadTop1Constit)->size() > 0)
	      std::cout << (*hadTop1Constit)->at(0).Pt() << "\t" << (*hadTop1Constit)->at(1).Pt() << "\t" << (*hadTop1Constit)->at(2).Pt() << std::endl;
	    if((*hadTop2Constit)->size() > 0)
	    std::cout << (*hadTop2Constit)->at(0).Pt() << "\t" << (*hadTop2Constit)->at(1).Pt() << "\t" << (*hadTop2Constit)->at(2).Pt() << std::endl;
	    if((*hadTop3Constit)->size() > 0)
	    std::cout << (*hadTop3Constit)->at(0).Pt() << "\t" << (*hadTop3Constit)->at(1).Pt() << "\t" << (*hadTop3Constit)->at(2).Pt() << std::endl;
	    HOTEval.printMatrixR(1);
	    HOTEval.printMatrixR(2);
	    HOTEval.printMatrixR(3);
	    HOTEval.printDimensions();
	    HOTEval.evaluateR();
	    HOTEval.printMatrixR(1);
	    HOTEval.printMatrixR(2);
	    HOTEval.printMatrixR(3);
	    for(const uint topflag: **FlagTop){
	      if(debug1) printf("\n\ttopflag = %6d", topflag);
	      if(topflag < 9999){
		if(topflag > 1020){
		  RecoTypes->Fill(2);
		  nRecoTops += 1;
		  flagcounter ++;
		  switch (flagcounter) {
		  case 1: Top1Reco = true;
		    break;
		  case 2: Top2Reco = true;
		    break;
		  case 3: Top3Reco = true;
		    break;
		  }
		if(debug1) printf("Booleans for Top reconstruction: %2d %2d %2d %2d", flagcounter, Top1Reco, Top2Reco, Top3Reco);
		if(debug1) printf("\n\t\tCnt: %2d \t Top1Reco: %2d \t Top2Reco: %2d \t Top3Reco: %2d", flagcounter, Top1Reco, Top2Reco, Top3Reco);
		}
		else if(topflag == 876 || topflag == 1004 || topflag == 1020){
		  RecoTypes->Fill(1);
		  nRecoWs += 1; //b not reconstructible, but q1 and q2 are
		}
		else
		  RecoTypes->Fill(0);
	      }
	    }
	    printf("\n\t\t\t# Fully Reconstuctible Tops: %2d \t Reconstructible W: %2d\n", nRecoTops, nRecoWs);

	    

            //print top properties
            for(const TopObject* top : tops)
            {
                //print basic top properties (top->p() gives a TLorentzVector)
                //N constituents refers to the number of jets included in the top
                //3 for resolved tops 
                //2 for W+jet tops
                //1 for fully merged AK8 tops
	      if(!silentRunning) printf("\tTop properties: Type: %3d,   Pt: %6.1lf,   Eta: %7.3lf,   Phi: %7.3lf,   M: %7.3lf,   Disc: %7.3f\n", static_cast<int>(top->getType()), top->p().Pt(), top->p().Eta(), top->p().Phi(), top->p().M(), top->getDiscriminator());

                //get vector of top constituents 
                const std::vector<Constituent const *>& constituents = top->getConstituents();

		//counter for jet matching
		uint mat1 = 0;
		uint mat2 = 0;
		uint mat3 = 0;
		uint mat4 = 0;

		//Matrix for jet matching
		uint tMatrixHOT[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
		

                //Print properties of individual top constituent jets 
                for(const Constituent* constituent : constituents)
                {
                    if(!silentRunning) printf("\t\tConstituent properties: Constituent type: %3d,   Pt: %6.1lf,   Eta: %7.3lf,   Phi: %7.3lf\n", constituent->getType(), constituent->p().Pt(), constituent->p().Eta(), constituent->p().Phi());

		    //Here, we are become gods, matching TLorentzVectors 
		    uint yui = 0;
		    uint yu2 = 0;
		    uint yu3 = 0;
		    uint yu4 = 0;
		    //testing that ->p() gives TLV that can exactly match original jet collection
		    // for(const TLorentzVector inJet : **AK4JetLV){
		    //   printf("\t\t\tJet %3d, Pt: %6.1lf, Match: %3d\n", yui++, inJet.Pt(), (inJet  == constituent->p()) );
		    //   mat1 += (inJet == constituent->p());
		    // }

		    //printf("\n");
		    for(const TLorentzVector inJet : **hadTop1Constit){
		      //printf("\t\t\thT1 %3d, Pt: %6.1lf, Match: %3d\n", yu2++, inJet.Pt(), (inJet  == constituent->p()) );
		      mat2 += (inJet == constituent->p());
		    }
		    //printf("\n");
		    for(const TLorentzVector inJet : **hadTop2Constit){
		      //printf("\t\t\thT2 %3d, Pt: %6.1lf, Match: %3d\n", yu3++, inJet.Pt(), (inJet  == constituent->p()) );
		      mat3 += (inJet == constituent->p());
		    }
		    //printf("\n");
		    for(const TLorentzVector inJet : **hadTop3Constit){
		      //printf("\t\t\thT3 %3d, Pt: %6.1lf, Match: %3d\n", yu4++, inJet.Pt(), (inJet  == constituent->p()) );
		      mat4 += (inJet == constituent->p());
		    }
                }    
		printf("\t\t\ttop1: %2d || top2: %2d|| top3: %2d || Discriminant: %6.4lf\n", mat2, mat3, mat4, top->getDiscriminator());
		//if(mat2 == 3 || mat3 == 3 || mat4 == 3) h_typeIII->Fill(top->getDiscriminator());
		//if(mat2 == 2 || mat3 == 2 || mat4 == 2) TypeWDiscr->Fill(top->getDiscriminator());
		//if(mat2 < 2 && mat3 < 2 && mat4 < 2) h_typeI->Fill(top->getDiscriminator());
            }

            //Print properties of the remaining system
            //the remaining system is used as the second portion of the visible system to calculate MT2 in the NT = 1 bin
            //const TopObject& rsys = ttr.getRsys();
            //printf("\tRsys properties: N constituents: %3d,   Pt: %6.1lf,   Eta: %7.3lf,   Phi: %7.3lf\n", rsys.getNConstituents(), rsys.p().Pt(), rsys.p().Eta(), rsys.p().Phi());
        
            printf("\n");
        }
    }
    catch(const TTException& e)
    {
        //Print exception message
        e.print();
        printf("Terminating run\n");
        fflush(stdout);

        exit(1);
    }

    //Open output file for histograms and tuples
    of = new TFile("results.root", "RECREATE");
    //Write HOT histos
    h_typeIII_hot->Write();
    h_typeIIb_hot->Write();
    h_typeIIw_hot->Write();
    h_typeIIo_hot->Write();
    h_typeI_hot->Write();
    h_type0_hot->Write();

    h_eventN1_III_hot->Write();
    h_eventN1_II_hot->Write();
    h_eventN1_I_hot->Write();
    h_eventN1_0_hot->Write();

    h_eventN2_III_hot->Write();
    h_eventN2_II_hot->Write();
    h_eventN2_I_hot->Write();
    h_eventN2_0_hot->Write();

    h_eventN3_III_hot->Write();
    h_eventN3_II_hot->Write();
    h_eventN3_I_hot->Write();
    h_eventN3_0_hot->Write();

    //Write BDT histos
    h_typeIII_bdt->Write();
    h_typeIIb_bdt->Write();
    h_typeIIw_bdt->Write();
    h_typeIIo_bdt->Write();
    h_typeI_bdt->Write();
    h_type0_bdt->Write();

    h_eventN1_III_bdt->Write();
    h_eventN1_II_bdt->Write();
    h_eventN1_I_bdt->Write();
    h_eventN1_0_bdt->Write();

    h_eventN2_III_bdt->Write();
    h_eventN2_II_bdt->Write();
    h_eventN2_I_bdt->Write();
    h_eventN2_0_bdt->Write();

    h_eventN3_III_bdt->Write();
    h_eventN3_II_bdt->Write();
    h_eventN3_I_bdt->Write();
    h_eventN3_0_bdt->Write();

    //Write others, write file, close
    // h_typeIII->Write();
    // h_typeIIb->Write();
    // h_typeIIw->Write();
    // h_typeIIo->Write();
    // TypeWDiscr->Write();
    // h_typeI->Write();

    RecoTypes->Write();
    of->Write();
    of->Close();
    
    //clean up pointers 
    delete AK4JetLV;
    delete AK4JetBtag;
    delete AK4qgMult;
    delete AK4qgPtD;
    delete AK4qgAxis1;
    delete AK4qgAxis2;
    delete AK4recoJetschargedHadronEnergyFraction;
    delete AK4recoJetschargedEmEnergyFraction;
    delete AK4recoJetsneutralEmEnergyFraction;
    delete AK4ElectronEnergyFraction;
    delete AK4PhotonEnergyFraction;
    delete AK4recoJetsneutralEnergyFraction;
    delete AK4recoJetsHFHadronEnergyFraction;
    delete AK4recoJetsmuonEnergyFraction;
    delete AK4recoJetsHFEMEnergyFraction;
    delete AK4NeutralHadronMultiplicity;
    delete AK4ChargedHadronMultiplicity;
    delete AK4ElectronMultiplicity;
    delete AK4MuonMultiplicity;
    delete AK4PhotonMultiplicity;
    delete AK4DeepCSVbb;
    delete AK4DeepCSVb;
    delete AK4DeepCSVc;
    delete AK4DeepCSVcc;
    delete AK4DeepCSVl;

    // delete AK8JetLV;
    // delete AK8SubjetLV;
    // delete AK8JetDeepAK8Top;
    // delete AK8JetSoftdropMass;

    //Delete histograms
    delete h_typeIII_hot;
    delete h_typeIIb_hot;
    delete h_typeIIw_hot;
    delete h_typeIIo_hot;
    delete h_typeI_hot;
    delete h_type0_hot;

    delete h_eventN1_III_hot;
    delete h_eventN1_II_hot;
    delete h_eventN1_I_hot;
    delete h_eventN1_0_hot;

    delete h_eventN2_III_hot;
    delete h_eventN2_II_hot;
    delete h_eventN2_I_hot;
    delete h_eventN2_0_hot;

    delete h_eventN3_III_hot;
    delete h_eventN3_II_hot;
    delete h_eventN3_I_hot;
    delete h_eventN3_0_hot;

    //Write BDT histos
    delete h_typeIII_bdt;
    delete h_typeIIb_bdt;
    delete h_typeIIw_bdt;
    delete h_typeIIo_bdt;
    delete h_typeI_bdt;
    delete h_type0_bdt;

    delete h_eventN1_III_bdt;
    delete h_eventN1_II_bdt;
    delete h_eventN1_I_bdt;
    delete h_eventN1_0_bdt;

    delete h_eventN2_III_bdt;
    delete h_eventN2_II_bdt;
    delete h_eventN2_I_bdt;
    delete h_eventN2_0_bdt;

    delete h_eventN3_III_bdt;
    delete h_eventN3_II_bdt;
    delete h_eventN3_I_bdt;
    delete h_eventN3_0_bdt;

    delete tf, tf2, td, tree, of;
    //delete h_typeIII, h_typeIIb, h_typeIIw, h_typeIIo, h_typeI, RecoTypes;
    exit(0);
}
