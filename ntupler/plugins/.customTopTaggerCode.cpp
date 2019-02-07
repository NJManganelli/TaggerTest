#include <cstdio>
#include <vector>

//File and plotting includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"

//manditory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
//this include is useful to get the helper function to make the vector of constituents
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"

//this include is necessary to handle exceptions thrown by the top tagger code
#include "TopTagger/CfgParser/include/TTException.h"

#include "rootdict.h"

//Includes to use BDT (based on TMVA)
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

//Struct for holding BDT candidate information
struct BDTCand {
  //Jet Indices
  uint idx_i;
  uint idx_j;
  uint idx_k;
  //BDT Variables, 2016 default...
  float btag;
  float ThPtOverSumPt;
  float AngleThWh;
  float AngleThBh;
  float HadrWmass;
  float TopMass;
  //BDT discriminant
  double discriminant;
};

double scaleBDT(double input){
  //return 1.666666666666*(input + 0.4);
  return 1.85*(input + 0.34);
}
double unscaleBDT(double input){
  //return 1.666666666666*(input + 0.4);
  return (input/1.85 - 0.34);
}

//Class for doing jet combinatorics and returning the list of best resolved top-tagged jets using the BDT method from the 2016 Four-top DL/SL analysis
class ResTTPermuter{
public:
  ResTTPermuter(std::vector<TLorentzVector>* inJets, std::vector<float, std::allocator<float> >*& inBtag);
  ResTTPermuter(std::vector<TLorentzVector>* inJets, std::vector<float, std::allocator<float> >*& inBtag, bool debug);
  ~ResTTPermuter();
  //uint GetNumPermutations();
  uint GetNumTrijets();
  bool ShouldSkip(uint jetIndex);
  std::vector<BDTCand> PitchCandidates(int candIndex);
  void CatchCandidates(std::vector<BDTCand> tossedCandidates); //should evaluate and store the best in the finalCandidates, plus store jet indices to be skipped next round
  void CatchCandidate(BDTCand tossedCandidate);
  std::vector<BDTCand> GetFinalIDXCandidates();
  std::vector<std::pair<std::vector<TLorentzVector>, double> > GetFinalCandidates();

private:
  bool _debug = false;
  uint _numTrijets = 0;
  uint _numPermutations = 0;
  uint _size = 0;
  //index list of jets that belong to a best candidate
  std::vector<uint> _skipIndices;
  //pointer to the AK4 jets and CSVv2 btags to be evaluated with the BDT method
  std::vector<TLorentzVector>* _inJets;
  std::vector<float> _inBtags;
  //Permit limitless number of candidates
  std::vector<std::vector<BDTCand> > _candVectorVector;
  //Create final candidates, storing the jet indices and the highest BDT discriminant from the 3 permutations (3 unique W-candidate jet-pairs) for the final best candidate
  std::vector<BDTCand> _finalCandidates;
};
//ResTTPermuter::ResTTPermuter(std::vector<TLorentzVector>* inJets, std::vector<double> inBtags){
ResTTPermuter::ResTTPermuter(std::vector<TLorentzVector>* inJets, std::vector<float, std::allocator<float> >*& inBtags){
  _inJets = inJets;
  _size = _inJets->size();
  //_inBtags = &*inBtags->begin();//Doesn't work... just iterate and assign
  for(std::vector<float>::iterator b = inBtags->begin(); b != inBtags->end(); b++)
    _inBtags.push_back(*b);
  int availJets = _inJets->size() - _skipIndices.size();
  _numTrijets = floor(availJets/3);
  // if(_numTrijets > 0)
  //   _canPermute = true;

}
ResTTPermuter::ResTTPermuter(std::vector<TLorentzVector>* inJets, std::vector<float, std::allocator<float> >*& inBtags, bool debug){
  _debug = debug;
  _inJets = inJets;
  _size = _inJets->size();
  for(std::vector<float>::iterator b = inBtags->begin(); b != inBtags->end(); b++)
    _inBtags.push_back(*b);
  int availJets = _inJets->size() - _skipIndices.size();
  _numTrijets = floor(availJets/3);
}
ResTTPermuter::~ResTTPermuter(){
}
uint ResTTPermuter::GetNumTrijets(){
  return _numTrijets;
}
// uint ResTTPermuter::GetNumPermutations(){
//   inline uint Factorial(uint in){
//     if(in < 2)
//       return 1;
//     else
//       return in*Factorial((in-1));
//   }
//   int availJets = (_inJets->size() - _skipIndices.size());
//   if(availJets > 2)
    

//     return 0;
//   else{
//     _canPermute = true;
//     return floor(possibilities/3);
//   }
// }
//Pitch all candidates, to be analyzed by BDT tagger, then Caught with updated discriminant and the best chosen
bool ResTTPermuter::ShouldSkip(uint jetIndex){
  for(int check = 0; check < _skipIndices.size(); check++)
    if(jetIndex == _skipIndices[check])
      return true; //if passed jet index matches any in the skip list
  return false; //if jet index doesn't match against any in the skip list
}
std::vector<BDTCand> ResTTPermuter::PitchCandidates(int candIndex){
  if(_debug)
    std::cout << "Starting to Pitch Candidates" << std::endl;
  std::vector<BDTCand> tempVector;
  if(candIndex < _numTrijets){
    for(uint i = 0; i < _size; i++){
      if(this->ShouldSkip(i))
	continue;
      for(uint j = i+1; j < _size; j++){
	if(this->ShouldSkip(j))
	  continue;
	for(uint k = j+1; k < _size; k++){
	  if(this->ShouldSkip(k))
	    continue;
	  TLorentzVector Wh, Bh, Th;
	  TLorentzVector first = (*_inJets)[i];
	  TLorentzVector second = (*_inJets)[j];
	  TLorentzVector third = (*_inJets)[k];

	  double angle12 = first.DeltaR(second);
	  double angle13 = first.DeltaR(third);
	  double angle23 = second.DeltaR(third);
	  double minAngle = fmin(fmin(angle12, angle13), angle23);
	  double theBtag = -9876.54321;
	  //Define the W hadron as the jet pair with minimum DeltaR
	  if(minAngle == angle12){
	    Wh = first + second;
	    Bh = third;
	    theBtag = _inBtags[k];
	  }
	  else if(minAngle == angle13){
	    Wh = first + third;
	    Bh = second;
	    theBtag = _inBtags[j];
	  }
	  else if(minAngle == angle23){
	    Wh = second + third;
	    Bh = first;
	    theBtag = _inBtags[i];
	  }
	  else
	    throw std::runtime_error("Logic Error! ResTTPermuter::PitchCandidates tripped on a Banana Peel");
	  Th = Wh + Bh;

	  BDTCand temp;
	  temp.idx_i = i;
	  temp.idx_j = j;
	  temp.idx_k = k;
	  temp.btag = theBtag;
	  temp.ThPtOverSumPt = Th.Pt()/(first.Pt() + second.Pt() + third.Pt());
	  temp.AngleThWh = fabs(Th.DeltaPhi(Wh)); //DeltaPhi separation
	  temp.AngleThBh = fabs(Th.DeltaPhi(Bh)); //DeltaPhi separation
	  temp.HadrWmass = Wh.M();
	  temp.TopMass = Th.M();
	  temp.discriminant = -9876.54321; //To be filled by TMVA and batted back to ResTTPermuter
	  tempVector.push_back(temp);
	}
      }
    }
  }
  if(_debug)
    std::cout << "\nSize of tempVector is ... " << tempVector.size() << std::endl;
  //_candVectorVector.push_back(tempVector);
  return tempVector;
}
void ResTTPermuter::CatchCandidates(std::vector<BDTCand> tossedCandidates){
  _finalCandidates.push_back(tossedCandidates[0]);
  _skipIndices.push_back(tossedCandidates[0].idx_i);
  _skipIndices.push_back(tossedCandidates[0].idx_j);
  _skipIndices.push_back(tossedCandidates[0].idx_k);
}
void ResTTPermuter::CatchCandidate(BDTCand tossedCandidate){
  _finalCandidates.push_back(tossedCandidate);
  _skipIndices.push_back(tossedCandidate.idx_i);
  _skipIndices.push_back(tossedCandidate.idx_j);
  _skipIndices.push_back(tossedCandidate.idx_k);
}
std::vector<BDTCand> ResTTPermuter::GetFinalIDXCandidates(){
  return _finalCandidates;
}
std::vector<std::pair<std::vector<TLorentzVector>, double> > ResTTPermuter::GetFinalCandidates(){
  std::vector<std::pair<std::vector<TLorentzVector>, double> > tempVector;
  for(uint finals = 0; finals < _numTrijets; finals++){
    std::pair<std::vector<TLorentzVector>, double> tempPair;
    std::vector<TLorentzVector> tempTLV;
    tempTLV.push_back( (*_inJets)[ _finalCandidates[finals].idx_i] );
    tempTLV.push_back( (*_inJets)[ _finalCandidates[finals].idx_j] );
    tempTLV.push_back( (*_inJets)[ _finalCandidates[finals].idx_k] );
    tempPair.first = tempTLV;
    tempPair.second = _finalCandidates[finals].discriminant;
    tempVector.push_back(tempPair);
    if(_debug)
      std::cout << "Chosen Indices: " << _finalCandidates[finals].idx_i << " " << _finalCandidates[finals].idx_j << " " << _finalCandidates[finals].idx_k << std::endl;
  }
  return tempVector;
}

class ResTTEvaluator{
 public:
  ResTTEvaluator(std::string topTaggerName);
  ResTTEvaluator(std::string topTaggerName, bool verbose, bool debug);
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
  void addReco(std::vector<TLorentzVector>* reco, std::vector<int> flagVector); //using boost::(type)_bitset would be more space efficient ... up to factor of 8
  void addReco(std::vector<TLorentzVector>* reco);
  std::pair< std::vector<TLorentzVector>, std::vector<int> > getReco(int index);
  void addAllReco( std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > *allRecoInput);

  void addGen(std::vector<TLorentzVector>* gen, std::vector<int> flagVector); //using boost::(type)_bitset would be more space efficient ... up to factor of 8
  void addGen(std::vector<TLorentzVector>* gen);
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
  
  //These methods do the classifying. They are called by evaluate as a submethod...
  void classifyR();
  void classifyG();

  //method for getting final pairs of 
  std::pair<double, std::string> getClassR(int index);
  std::vector<std::pair<double, std::string>> getClassesR();
  std::vector<std::pair<double, std::string>> getOrderedClassesR();

  std::vector<std::pair<double, int>> getOrderedAggregatesR();

 private:
  bool _debug = false;
  bool _verbose = false;
  const static uint _maxSize = 100;
  uint _nReco;
  uint _nGen;
  uint _nCand;
  uint _matchR[_maxSize][_maxSize][3] = {0}; //support max 100 candidates and gen tops... shouldn't ever be a problem at 13TeV CoM Energy 
  uint _matchG[_maxSize][_maxSize][3] = {0}; //support max 100 candidates and gen tops... shouldn't ever be a problem at 13TeV CoM Energy 
  int bestHMatchR[_maxSize];
  int bestNumBR[_maxSize];
  int bestNumQR[_maxSize];
  int totNumBR[_maxSize];
  int totNumQR[_maxSize];

  //encode in this matrix the info for whether it was even POSSIBLE to match that jet... 
  bool _haveEvaluatedR;
  bool _haveClassifiedR;
  bool _haveEvaluatedG;
  bool _haveClassifiedG;
  std::string _topTaggerName;
  uint _orderIndex[_maxSize];
  std::vector< std::pair< std::vector<TLorentzVector>, double > > _cand;
  std::vector< std::pair< std::vector<TLorentzVector>, double >* > _ordCand; //use pointers?
  std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > _reco; //include flags in here...
  std::vector< std::pair< std::vector<TLorentzVector>, std::vector<int> > > _gen; //include flags in here...
  std::vector< std::pair< double, std::string > > _classR; //classification, using number and string
  std::vector< std::pair< double, std::string > > _ordClassR; //shoulda made a struct...
  std::vector< std::pair< double, int > > _aggrJetsR; //store total number of matched jets
  std::vector< std::pair< double, int > > _ordAggrJetsR; //store in ordered fashion, to be retrieved simultaneously...
  bool _haveFlags;
};
ResTTEvaluator::ResTTEvaluator(std::string topTaggerName){
  _nReco = 0;
  _nGen = 0;
  _nCand = 0;
  _haveEvaluatedR = false;
  _haveClassifiedR = false;
  _haveEvaluatedG = false;
  _haveClassifiedG = false;
  _haveFlags = false;
  _topTaggerName = topTaggerName;
  //std::cout << "ResTTEvaluator \"" << _topTaggerName << "\" being created" << std::endl;
}
ResTTEvaluator::ResTTEvaluator(std::string topTaggerName, bool verbose, bool debug){
  _verbose = verbose;
  _debug = debug;
  _nReco = 0;
  _nGen = 0;
  _nCand = 0;
  _haveEvaluatedR = false;
  _haveClassifiedR = false;
  _haveEvaluatedG = false;
  _haveClassifiedG = false;
  _haveFlags = false;
  _topTaggerName = topTaggerName;
  if(_verbose)
    std::cout << "ResTTEvaluator \"" << _topTaggerName << "\" being created" << std::endl;
}
ResTTEvaluator::~ResTTEvaluator(){
  if(_verbose)
    std::cout << "ResTTEvaluator \"" << _topTaggerName << "\" being destroyed" << std::endl;
}
void ResTTEvaluator::addCand(std::vector<TLorentzVector>* cand, double discriminant){
  std::pair< std::vector<TLorentzVector>, double> tempCand;
  tempCand.first = *cand;
  tempCand.second = discriminant;
  if(_verbose){
    std::cout << "\nAdding Cand: " << std::endl;
    for(int q = 0; q < 3; q++)
      std::cout << "\t" << (tempCand.first[q]).Pt() << "\t" << (tempCand.first[q]).Eta() 
	      << "\t" << (tempCand.first[q]).Phi() << "\t" << (tempCand.first[q]).E() << std::endl;
  }
  _cand.push_back(tempCand);
  _nCand++;
}
std::vector<std::pair<std::vector<TLorentzVector>, double>> ResTTEvaluator::getAllCand(){
  return _cand;
}
std::pair< std::vector<TLorentzVector>, double > ResTTEvaluator::getCand(int index){
  if(index < _nCand)
    return _cand[index];
  else
    std::cout << "This index is not in range" << std::endl;
}
double ResTTEvaluator::getDiscr(int index){
  return (_cand[index]).second;
}
std::pair< std::vector<TLorentzVector>, double > ResTTEvaluator::getOrderedCand(int index){
  std::cout << "Method not implemented, just returning unordered candidate..." << std::endl;
  if(index < _nCand)
    return _cand[index];
  else
    std::cout << "This index is not in range" << std::endl;
}
double ResTTEvaluator::getOrderedDiscr(int index){
  std::cout << "Method not implemented, just returning unordered candidate's discriminant..." << std::endl;
  if(index < _nCand)
    return (_cand[index]).second;
  else
    std::cout << "This index is not in range" << std::endl;
}
void ResTTEvaluator::addReco(std::vector<TLorentzVector>* reco, std::vector<int> flags){
  if(_reco.size() > 0 && !_haveFlags){
    _haveFlags = false; //explicit for readibility
    std::cout << "Added genReco with flags, but previous candidate didn't have flags! Will not perform any calculations depending on presence of flags..." << std::endl;
  }
  else
    _haveFlags = true;
  std::pair< std::vector<TLorentzVector>, std::vector<int> > tempReco;
  tempReco.first = *reco;
  tempReco.second = flags;
  if(_verbose){
    std::cout << "\nAdding Reco: ";
    for(int qq = 0; qq < flags.size(); qq++)
      std::cout << "\t\tfl" << qq << ": " << flags[qq];
    std::cout << std::endl;
    for(int q = 0; q < 3; q++)
      std::cout << "\t" << (tempReco.first[q]).Pt() << "\t" << (tempReco.first[q]).Eta() 
		<< "\t" << (tempReco.first[q]).Phi() << "\t" << (tempReco.first[q]).E() << std::endl;
  }
  _reco.push_back(tempReco);
  _nReco++;
}
void ResTTEvaluator::addReco(std::vector<TLorentzVector>* reco){
  if(_haveFlags){
    _haveFlags = false;
    std::cout << "Added genReco without flags, but previous candidate had flags! Will not perform any calculations depending on presence of flags..." << std::endl;
  }
  std::vector<int> tempFlags;
  tempFlags.push_back(-1);
  std::pair<std::vector<TLorentzVector>, std::vector<int> > tempReco;
  tempReco.first = *reco;
  tempReco.second = tempFlags;
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
  if(index < _nReco)
    return _reco[index];
  else
    std::cout << "\n\tThis index is not in range!" << std::endl;
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
  if(_verbose){
    std::cout << "\nAdding Gen: " << std::endl;
  for(int q = 0; q < 3; q++)
    std::cout << "\t" << (tempGen.first[q]).Pt() << "\t" << (tempGen.first[q]).Eta() 
		<< "\t" << (tempGen.first[q]).Phi() << "\t" << (tempGen.first[q]).E() << std::endl;
  }
  _gen.push_back(tempGen);
  _nGen++;
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
  if(index < _nGen)
    return _gen[index];
  else
    std::cout << "\n\tThis index is not in range!" << std::endl;
}
void ResTTEvaluator::printCand(int index){
  std::cout << "This isn't implemented yet... " << std::endl;
}
void ResTTEvaluator::printGen(int index){
  std::cout << "This isn't implemented yet... " << std::endl;
}
void ResTTEvaluator::printMatrixR(int index){
  //std::cout << "index: " << index << " nReco: " << _nReco << " (-1 < index && index < nReco): " << (-1 < index && index < _nReco) << std::endl;
  if(-1 < index && index < _nReco){ //could also do negative indices for searching from the end... hmmm
    if(!_haveEvaluatedR)
      std::cout << "\nPre-evaluation matrix" << std::endl;
    std::cout << "\n\t Match Matrix for Top Candidate " << index << "\n\t\t\tb jet\tq1 jet\tq2 jet";
    for(int i = 0; i < _nReco; i++)
      std::cout << "\n\t Gen Reco " << i << "\t" << _matchR[index][i][0] << "\t" 
		<< _matchR[index][i][1] << "\t" << _matchR[index][i][2];
    std::cout << std::endl;
  }
  else
    std::cout << "\n\tThis index is not in range!" << std::endl;
}
void ResTTEvaluator::printMatrixG(int index){
  if(-1 < index && index < _nGen){
    if(!_haveEvaluatedG)
      std::cout << "\nPre-evaluation matrix" << std::endl;
    std::cout << "\n\t Match Matrix for Top Candidate " << index << "\n\t\t\tb jet\tq1 jet\tq2 jet";
    for(int i = 0; i < _nGen; i++)
      std::cout << "\n\t Gen Part " << i+1 << "\t" << _matchG[index][i][0] << "\t" 
		<< _matchG[index][i][1] << "\t" << _matchG[index][i][2];
    std::cout << std::endl;
  }
  else
    std::cout << "\n\tThis index is not in range!" << std::endl;
}
void ResTTEvaluator::printDimensions(){
  std::cout << "nGen: " << _nGen << "\tnReco: " << _nReco << "\tnCand: " << _nCand << "\tmaxSize: " << _maxSize << std::endl;
}
void ResTTEvaluator::evaluateR(){
  for(int m = 0; m < _nCand; m++) //evaluate for each candidate
    for(int n = 0; n < _cand[m].first.size(); n++) //each candidate's constituent
      //for(int ii = 0; ii < 3; ii++) //each candidate's constituent
      for(int o = 0; o < _nReco; o++) //for each candidate, check each Reconstructed top quark's constituents
	for(int p = 0; p < _reco[o].first.size(); p++) //more permissive...
	//for(int jj = 0; jj < 3; jj++) //check all the jets in the collection... should be 3 at all times... need to protect...
	  {
	    if(_debug)
	      std::cout << "\n\ti|ii|j|jj: " << m << n << o << p << "\t" << _cand[m].first[n].Pt() << "\t" << _reco[o].first[p].Pt() << "\t" << m << o << p;
	    int test = static_cast<int>(_reco[o].first[p] == _cand[m].first[n]);
	    _matchR[m][o][p] = (test > _matchR[m][o][p] ? test : _matchR[m][o][p]); //for each candidate, store match truth in the matrix
	    if(_debug)
	      std::cout << "\t" << test << "\t" << _matchR[m][o][p];
	  }
  this->classifyR();
  _haveEvaluatedR = true;
  std::cout << std::endl;
}
void ResTTEvaluator::evaluateG(){
  std::cout << "If I were a real little evaluator (Gen), I would have done some evaluation (and stuff!). Since I'm not (yet), I'll just let you know this worked!" << std::endl;
}
void ResTTEvaluator::classifyR(){

  int HMatch[_nCand][_nReco] = {0};
  int sumHMatch[_nCand] = {0};
  // int VMatch_b[_nCand] = {0};
  // int VMatch_q1[_nCand] = {0};
  // int VMatch_q2[_nCand] = {0};
  for(int m = 0; m < _nCand; m++){
    bestHMatchR[m] = 0;
    bestNumBR[m] = 0;
    bestNumQR[m] = 0;
    totNumBR[m] = 0;
    totNumQR[m] = 0;
    //Horizontal matches
    for(int o = 0; o < _nReco; o++){
      int thisHMatch = 0;
      for(int p = 0; p < 3; p++)
	thisHMatch += _matchR[m][o][p];
      if(_debug)
	std::cout << "\n\t\t\t\tmo: " << m << o << " best: " << bestHMatchR[m] << " this: " << thisHMatch;
      bool bettermatch = (thisHMatch > bestHMatchR[m]);
      bestHMatchR[m] = (bettermatch ? thisHMatch : bestHMatchR[m]);
      if(_debug)
	std::cout << " best: " << bestHMatchR[m];
      bestNumBR[m] = (bettermatch ? _matchR[m][o][0] : bestNumBR[m]);
      bestNumQR[m] = (bettermatch ? (_matchR[m][o][1] + _matchR[m][o][2]) : bestNumQR[m]);
      // bestHMatchR[m] = (thisHMatch > bestHMatchR[m] ? thisHMatch : bestHMatchR[m]);
      // bestNumBR[m] = (thisHMatch > bestHMatchR[m] ? _matchR[m][o][0] : bestNumBR[m]);
      // bestNumQR[m] = (thisHMatch > bestHMatchR[m] ? (_matchR[m][o][1] + _matchR[m][o][2]) : bestNumQR[m]);
      sumHMatch[m] += thisHMatch;
    }
    //Vertical Match
    //for(int p = 0; p < 3; p++){
    // int thisVMatch_b = 0;
    // int thisVMatch_q1 = 0;
    // int thisVMatch_q2 = 0;
    for(int o = 0; o < _nReco; o++){
      totNumBR[m] += _matchR[m][o][0];
      totNumQR[m] += _matchR[m][o][1];
      totNumQR[m] += _matchR[m][o][2];
    }
    if(_verbose)
      std::cout << "\n\t\t\t" << m << "\tBestB: " << bestNumBR[m] << "\tBestQ: " << bestNumQR[m] 
		<< "\tTotB: " << totNumBR[m] << "\tTotQ: " << totNumQR[m] << std::endl;

    //UNHOLY WORK
    //final classification
    //Test for constraint failure:
    if(bestNumBR[m] + bestNumQR[m] > 3 || totNumBR[m] + totNumQR[m] > 3)
      std::cout << "WARNING: Constraint failure, best or total number of jet matches exceeds 3" << std::endl;
    std::string tempClass = "";
    //best case had 2 q jets matched
    if(bestNumQR[m] == 2){
      //best case had 1 b jet matched
      if(bestNumBR[m] == 1){
	//perfect match
	tempClass = "typeIII";
      }
      //best case had 0 b jet matched
      else if(bestNumBR[m] == 0){
	//there was 1 b jet matched
	if(totNumBR[m] == 1){
	  //swapped a b with another Reco
	  tempClass = "typeIIb";
	}
	//there were 0 b jets matched
	else if(totNumBR[m] == 0){
	  //there were 3 q jets matched
	  if(totNumQR[m] == 3){
	    //misidentified a q jet AS the 'b' in candidate
	    tempClass = "typeIImib";
	  }
	  //there were 2 q jets matched
	  else if(totNumQR[m] == 2){
	    //misidentified some non-reco/non-quark jet as the 'b' in candidate
	    tempClass = "typeIImib";
	  }
	  else
	    std::cout << "WARNING: LOGIC FAILURE!" << std::endl;
	}
	else
	  std::cout << "WARNING: LOGIC FAILURE!" << std::endl;
      }
      else
	std::cout << "WARNING: LOGIC FAILURE!" << std::endl;
    }
    //best case had 1 q jet matched
    else if(bestNumQR[m] == 1){
      //best case had 1 b jet matched
      if(bestNumBR[m] == 1){
	//there was 1 b jet matched
	if(totNumBR[m] == 1){
	  //there were 2 q jets matched
	  if(totNumQR[m] == 2){
	    //swapped q  with another Reco
	    tempClass = "typeIIw";
	  }
	  //there was 1 q jet matched
	  else if(totNumQR[m] == 1){
	    //misidentified some other jet AS the 'q' in candidate
	    tempClass = "typeIImiq";
	  }
	  else
	    //best # 1 = 1 and tot # q = 0 or more
	    std::cout << "WARNING: LOGIC FAILURE!" << std::endl;
	}
	//there were 2 b jets matched
	else if(totNumBR[m] == 2){
	  //there was 1 q jet matched
	  if(totNumQR[m] == 1){
	    //misidentify a b jet AS the 'q' in candidate
	    tempClass = "typeIImiq";
	  }
	  else
	    std::cout << "WARNING: LOGIC FAILURE!" << std::endl;
	}
	else
	  //best # b = 1 and tot # b = 0...
	  std::cout << "WARNING: LOGIC FAILURE!" << std::endl;

      }
      else if(bestNumBR[m] == 0){
	//identified only 1 q in the best Reco, no b in that Reco... could be unreconstructable jets?
	if( (totNumBR[m] + totNumQR[m]) == 3)
	  tempClass = "typeIt";
	else if( (totNumBR[m] + totNumQR[m]) == 2)
	  tempClass = "typeIp";
	else if( (totNumBR[m] + totNumQR[m]) == 1)
	  tempClass = "typeI";
	else
	  tempClass = "type0";
      }
      else
	std::cout << "WARNING: LOGIC FAILURE!" << std::endl;

    }
    //identified no q jets... amazingly!
    else if(bestNumQR[m] == 0){
      //found the b at least...
      if(bestNumBR[m] == 1){
	//found more than one b...
	if(totNumBR[m] > 1){
	  std::cout << "We may have found the triple 'b'entente!" << std::endl;
	  if( (totNumBR[m] + totNumQR[m]) == 3)
	    tempClass = "typeIt";
	  else if( (totNumBR[m] + totNumQR[m]) == 2)
	    tempClass = "typeIp";
	  else if( (totNumBR[m] + totNumQR[m]) == 1)
	    tempClass = "typeI";
	  else
	    tempClass = "type0";
	}
	else if(totNumBR[m] == 1){
	  if( (totNumBR[m] + totNumQR[m]) == 3)
	    tempClass = "typeIt";
	  else if( (totNumBR[m] + totNumQR[m]) == 2)
	    tempClass = "typeIp";
	  else if( (totNumBR[m] + totNumQR[m]) == 1)
	    tempClass = "typeI";
	  else
	    tempClass = "type0";
	}
	else
	  std::cout << "WARNING: LOGIC FAILURE!" << std::endl;
      }
      //no b jets, no q jets...
      else if(bestNumBR[m] == 0){
	if(_verbose)
	  std::cout << "We apparently found a perfectly UNMATCHED case... Type0. Complete Failure. Amazing!" << std::endl;
	tempClass = "type0";
      }
      else
	std::cout << "WARNING: LOGIC FAILURE!" << std::endl;
    }
    else
      std::cout << "WARNING: LOGIC FAILURE!" << std::endl;
    
    //store class
    std::pair< double, std::string> tempClassPair;
    tempClassPair.first = _cand[m].second;
    tempClassPair.second = tempClass;
    _classR.push_back(tempClassPair);
    std::pair< double, uint> tempAggrPair;
    tempAggrPair.first = _cand[m].second;
    tempAggrPair.second = totNumBR[m]+totNumQR[m];
    _aggrJetsR.push_back(tempAggrPair);
    if(_verbose)
      std::cout << "\t\t\t\tClassification: "  << _classR[m].first << "\t\t" << _classR[m].second << std::endl;
  }
  //for(int mm = 0; mm < _nCand; mm++)
    //std::cout << "\t\t\t\tClassification: "  << _classR[mm].first << "\t\t" << _classR[mm].second << std::endl;

  //make ordered classification vector
  _ordClassR = _classR;
  std::sort(_ordClassR.begin(), _ordClassR.end(), [](const std::pair<double, std::string> &left, const std::pair<double, std::string> &right) {
      return left.first > right.first;
    });
  _ordAggrJetsR = _aggrJetsR;
  std::sort(_ordAggrJetsR.begin(), _ordAggrJetsR.end(), [](const std::pair<double, int> &left, const std::pair<double, int> &right) {
      return left.first > right.first;
    });
  
}
void ResTTEvaluator::classifyG(){
  std::cout << "If I were a real little classifer (Gen), I would have done some classification (and stuff!). Since I'm not (yet), I'll just let you know this worked!" << std::endl;
}
std::pair<double, std::string> ResTTEvaluator::getClassR(int index){
  return _classR[index];
}
std::vector<std::pair<double, std::string>> ResTTEvaluator::getClassesR(){
  return _classR;
}
std::vector<std::pair<double, std::string>> ResTTEvaluator::getOrderedClassesR(){
  // std::sort(_ordClassR.begin(), _ordClassR.end(), [](const std::pair<;
  // 	    std::sort(v.begin(), v.end(), [](const std::pair<int,int> &left, const std::pair<int,int> &right) {
  // 		return left.second < right.second;
  // 	      });
  //could be moved to classify or evaluate method instead, so no double sorting
  // std::vector<std::pair<double, std::string>> _ordClassRout = _classR;
  // std::sort(_ordClassRout.begin(), _ordClassRout.end(), [](const std::pair<double, std::string> &left, const std::pair<double, std::string> &right) {
  //     return left.first > right.first;
  //   });
  // for(int rr = 0; rr < _ordClassR.size(); rr++){
  //   _ordClassRout.push_back(*(_ordClassR[rr]));
  // }
  //return _ordClassRout;
  return _ordClassR;
}
std::vector<std::pair<double, int>> ResTTEvaluator::getOrderedAggregatesR(){
  return _ordAggrJetsR;
}


int main(int argc, char* argv[])
{
  const std::string inFile = argv[1];
  const std::string outFile = argv[2];

  std::cout << "Input file: " << inFile << " Output file: " << outFile << std::endl;

  //bool for silencing original top quark properties (except event #)
  bool silentRunning = true;
  bool debug1 = false;
  bool debug2 = false;
  bool isNewStyle = true;
  //silentRunning = true;
  bool runStdExample = false;
  TFile *tf, *tf2, *of;
  TDirectory *td;
  TTree *tree;
  //std::string postfix = "tttt";
  uint maxEventsToProcess = -1;

  //nJet and candidate counting histos
  TH1I *h_nJet = new TH1I ("h_nJet", "Event Jet Distribution; nJets; Number of Events", 18, 6, 24);
  TH1I *h_nJet_N1_hot = new TH1I ("h_nJet_N1_hot", "Number of events with 1st HOT candidate; nJets; Number of Events", 18, 6, 24);
  TH1I *h_nJet_N2_hot = new TH1I ("h_nJet_N2_hot", "Number of events with 2nd HOT candidate; nJets; Number of Events", 18, 6, 24);
  TH1I *h_nJet_N3_hot = new TH1I ("h_nJet_N3_hot", "Number of events with 3rd HOT candidate; nJets; Number of Events", 18, 6, 24);
  TH1I *h_nJet_N1_bdt = new TH1I ("h_nJet_N1_bdt", "Number of events with 1st BDT candidate; nJets; Number of Events", 18, 6, 24);
  TH1I *h_nJet_N2_bdt = new TH1I ("h_nJet_N2_bdt", "Number of events with 2nd BDT candidate; nJets; Number of Events", 18, 6, 24);
  TH1I *h_nJet_N3_bdt = new TH1I ("h_nJet_N3_bdt", "Number of events with 3rd BDT candidate; nJets; Number of Events", 18, 6, 24);

  //Event Info Histograms
  TH1I *h_nTrueRecoTops_full = new TH1I ("h_nTrueRecoTops_full", 
					 "Number of true hadronic tops that are fully in event selection and acceptance; Tops per event; Events", 6, 0, 6);
  TH1I *h_nTrueRecoTops = new TH1I ("h_nTrueRecoTops", 
				    "Number of true hadronic tops in events surviving event selection; Tops per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_hot = new TH1I ("h_nCandTops_hot", 
				    "Number of tagger candidate (HOT) tops in events surviving event selection; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_bdt = new TH1I ("h_nCandTops_bdt", 
				    "Number of tagger candidate (BDT) tops in events surviving event selection; Top candidates per event; Events", 6, 0, 6);
  //All reco top counting
  TH1I *h_nCandTops_hot1 = new TH1I ("h_nCandTops_hot1", 
				     "Number of tagger candidate (HOT) tops in events with 1 real hadronic top; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_bdt1 = new TH1I ("h_nCandTops_bdt1", 
				     "Number of tagger candidate (BDT) tops in events with 1 real hadronic top; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_hot2 = new TH1I ("h_nCandTops_hot2", 
				     "Number of tagger candidate (HOT) tops in events with 2 real hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_bdt2 = new TH1I ("h_nCandTops_bdt2", 
				     "Number of tagger candidate (BDT) tops in events with 2 real hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_hot3 = new TH1I ("h_nCandTops_hot3", 
				     "Number of tagger candidate (HOT) tops in events with 3 real hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_bdt3 = new TH1I ("h_nCandTops_bdt3", 
				     "Number of tagger candidate (BDT) tops in events with 3 real hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_hot4 = new TH1I ("h_nCandTops_hot4", 
				     "Number of tagger candidate (HOT) tops in events with 4 real hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_bdt4 = new TH1I ("h_nCandTops_bdt4", 
				     "Number of tagger candidate (BDT) tops in events with 4 real hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_hot_oth = new TH1I ("h_nCandTops_hot_oth", 
				     "Number of tagger candidate (HOT) tops in events with other category; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandTops_bdt_oth = new TH1I ("h_nCandTops_bdt_oth", 
				     "Number of tagger candidate (BDT) tops in events with other category; Top candidates per event; Events", 6, 0, 6);
  //Full Reco top counting
  TH1I *h_nCandFullTops_hot_0 = new TH1I ("h_nCandFullTops_hot_0", 
				     "Number of tagger candidate (HOT) tops in events with 0 full-reco hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_bdt_0 = new TH1I ("h_nCandFullTops_bdt_0", 
				     "Number of tagger candidate (BDT) tops in events with 0 full-reco hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_hot1 = new TH1I ("h_nCandFullTops_hot1", 
				     "Number of tagger candidate (HOT) tops in events with 1 full-reco hadronic top; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_bdt1 = new TH1I ("h_nCandFullTops_bdt1", 
				     "Number of tagger candidate (BDT) tops in events with 1 full-reco hadronic top; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_hot2 = new TH1I ("h_nCandFullTops_hot2", 
				     "Number of tagger candidate (HOT) tops in events with 2 full-reco hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_bdt2 = new TH1I ("h_nCandFullTops_bdt2", 
				     "Number of tagger candidate (BDT) tops in events with 2 full-reco hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_hot3 = new TH1I ("h_nCandFullTops_hot3", 
				     "Number of tagger candidate (HOT) tops in events with 3 full-reco hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_bdt3 = new TH1I ("h_nCandFullTops_bdt3", 
				     "Number of tagger candidate (BDT) tops in events with 3 full-reco hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_hot4 = new TH1I ("h_nCandFullTops_hot4", 
				     "Number of tagger candidate (HOT) tops in events with 4 full-reco hadronic tops; Top candidates per event; Events", 6, 0, 6);
  TH1I *h_nCandFullTops_bdt4 = new TH1I ("h_nCandFullTops_bdt4", 
				     "Number of tagger candidate (BDT) tops in events with 4 full-reco hadronic tops; Top candidates per event; Events", 6, 0, 6);


  //HOT discriminant histos
  int nBins = 100;
  TH2F *h_typeIII_hot = new TH2F ("h_typeIII_hot_", "Type III (correct) Top Quarks; Discriminant; Number of Jets", 
				  nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIIb_hot = new TH2F ("h_typeIIb_hot_", "Type II (b swapped) Top Quarks; Discriminant; Number of Jets", 
				  nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIImib_hot = new TH2F ("h_typeIImib_hot_", "Type II (misidentified b from anywhere non-b) Top Quarks; Discriminant; Number of Jets", 
				    nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIIw_hot = new TH2F ("h_typeIIw_hot_", "Type II (q1 or q2 swapped) Top Quarks; Discriminant; Number of Jets", 
				  nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIImiq_hot = new TH2F ("h_typeIImiq_hot_", "Type II (misidentified q from anywhere non-q) Top Quarks; Discriminant; Number of Jets", 
				    nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIt_hot = new TH2F ("h_typeIt_hot_", "Type I (3 top-daughters matched, 1 per reco top); Discriminant; Number of Jets", 
				 nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_typeIp_hot = new TH2F ("h_typeIp_hot_", "Type I (2 top-daughters matched, 1 per reco top); Discriminant; Number of Jets", 
				 nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_type0xI_hot = new TH2F ("h_type0xI_hot_", "Type 0xI (all other tagger candidates); Discriminant ; Number of Jets", 
				  nBins, 0.0, 1.0, 7, 7, 14);

  TH2F *h_eventN1_sum_hot = new TH2F ("h_eventN1_sum_hot", "Highest Disc Cand (all) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_III_hot = new TH2F ("h_eventN1_tIII_hot", "Highest Disc Cand (Type III) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_II_hot = new TH2F ("h_eventN1_tII_hot", "Highest Disc Cand (Type II) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_IIb_hot = new TH2F ("h_eventN1_tIIb_hot", "Highest Disc Cand (Type IIb) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_IImib_hot = new TH2F ("h_eventN1_tIImib_hot", "Highest Disc Cand (Type IImib) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_IIw_hot = new TH2F ("h_eventN1_tIIw_hot", "Highest Disc Cand (Type IIw) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_IImiq_hot = new TH2F ("h_eventN1_tIImiq_hot", "Highest Disc Cand (Type IImiq) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_Itp_hot = new TH2F ("h_eventN1_tItp_hot", "Highest Disc Cand (Type Itp) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_It_hot = new TH2F ("h_eventN1_tIt_hot", "Highest Disc Cand (Type It) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_Ip_hot = new TH2F ("h_eventN1_tIp_hot", "Highest Disc Cand (Type Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_0xI_hot = new TH2F ("h_eventN1_t0x1_hot", "Highest Disc Cand (Type 0x1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_I_hot = new TH2F ("h_eventN1_t1_hot", "Highest Disc Cand (Type 1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_0_hot = new TH2F ("h_eventN1_t0_hot", "Highest Disc Cand (Type 0) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);

  TH2F *h_eventN2_sum_hot = new TH2F ("h_eventN2_sum_hot", "2nd Highest Disc Cand (all) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_III_hot = new TH2F ("h_eventN2_tIII_hot", "2nd Highest Disc Cand (Type III) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_II_hot = new TH2F ("h_eventN2_tII_hot", "2nd Highest Disc Cand (Type II) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_IIb_hot = new TH2F ("h_eventN2_tIIb_hot", "2nd Highest Disc Cand (Type IIb) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_IImib_hot = new TH2F ("h_eventN2_tIImib_hot", "2nd Highest Disc Cand (Type IImib) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_IIw_hot = new TH2F ("h_eventN2_tIIw_hot", "2nd Highest Disc Cand (Type IIw) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_IImiq_hot = new TH2F ("h_eventN2_tIImiq_hot", "2nd Highest Disc Cand (Type IImiq) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_Itp_hot = new TH2F ("h_eventN2_tItp_hot", "2nd Highest Disc Cand (Type Itp) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_It_hot = new TH2F ("h_eventN2_tIt_hot", "2nd Highest Disc Cand (Type It) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_Ip_hot = new TH2F ("h_eventN2_tIp_hot", "2nd Highest Disc Cand (Type Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_0xI_hot = new TH2F ("h_eventN2_t0x1_hot", "2nd Highest Disc Cand (Type 0x1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_I_hot = new TH2F ("h_eventN2_t1_hot", "2nd Highest Disc Cand (Type 1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_0_hot = new TH2F ("h_eventN2_t0_hot", "2nd Highest Disc Cand (Type 0) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);

  TH2F *h_eventN3_sum_hot = new TH2F ("h_eventN3_sum_hot", "3rd Highest Disc Cand (all) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_III_hot = new TH2F ("h_eventN3_tIII_hot", "3rd Highest Disc Cand (Type III) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_II_hot = new TH2F ("h_eventN3_tII_hot", "3rd Highest Disc Cand (Type II) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_IIb_hot = new TH2F ("h_eventN3_tIIb_hot", "3rd Highest Disc Cand (Type IIb) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_IImib_hot = new TH2F ("h_eventN3_tIImib_hot", "3rd Highest Disc Cand (Type IImib) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_IIw_hot = new TH2F ("h_eventN3_tIIw_hot", "3rd Highest Disc Cand (Type IIw) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_IImiq_hot = new TH2F ("h_eventN3_tIImiq_hot", "3rd Highest Disc Cand (Type IImiq) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_Itp_hot = new TH2F ("h_eventN3_tItp_hot", "3rd Highest Disc Cand (Type It / Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_It_hot = new TH2F ("h_eventN3_tIt_hot", "3rd Highest Disc Cand (Type It) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_Ip_hot = new TH2F ("h_eventN3_tIp_hot", "3rd Highest Disc Cand (Type Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_0xI_hot = new TH2F ("h_eventN3_t0x1_hot", "3rd Highest Disc Cand (Type 0 / I) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_I_hot = new TH2F ("h_eventN3_t1_hot", "3rd Highest Disc Cand (Type 1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_0_hot = new TH2F ("h_eventN3_t0_hot", "3rd Highest Disc Cand (Type 0) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);

  //BDT Candidates
  TH2F *h_typeIII_bdt = new TH2F ("h_typeIII_bdt_", "Type III (correct) Top Quarks; Discriminant; Number of Jets", 
				  nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIIb_bdt = new TH2F ("h_typeIIb_bdt_", "Type II (b swapped) Top Quarks; Discriminant; Number of Jets", 
				  nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIIw_bdt = new TH2F ("h_typeIIw_bdt_", "Type II (q1 or q2 swapped) Top Quarks; Discriminant; Number of Jets", 
				  nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIImib_bdt = new TH2F ("h_typeIImib_bdt_", "Type II (misidentified b from anywhere non-b) Top Quarks; Discriminant; Number of Jets", 
				    nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIImiq_bdt = new TH2F ("h_typeIImiq_bdt_", "Type II (misidentified q from anywhere non-q) Top Quarks; Discriminant; Number of Jets", 
				    nBins, 0.0, 1.0, 7, 7, 14); 
  TH2F *h_typeIt_bdt = new TH2F ("h_typeIt_bdt_", "Type I (3 top-daughters matched, 1 per reco top); Discriminant; Number of Jets", 
				 nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_typeIp_bdt = new TH2F ("h_typeIp_bdt_", "Type I (2 top-daughters matched, 1 per reco top); Discriminant; Number of Jets", 
				 nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_type0xI_bdt = new TH2F ("h_type0xI_bdt_", "Type 0x1 (all other tagger candidates); Discriminant ; Number of Jets", 
				  nBins, 0.0, 1.0, 7, 7, 14);

  TH2F *h_eventN1_sum_bdt = new TH2F ("h_eventN1_sum_bdt", "Highest Disc Cand (all) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_III_bdt = new TH2F ("h_eventN1_tIII_bdt", "Highest Disc Cand (Type III) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_II_bdt = new TH2F ("h_eventN1_tII_bdt", "Highest Disc Cand (Type II) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_IIb_bdt = new TH2F ("h_eventN1_tIIb_bdt", "Highest Disc Cand (Type IIb) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_IImib_bdt = new TH2F ("h_eventN1_tIImib_bdt", "Highest Disc Cand (Type IImib) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_IIw_bdt = new TH2F ("h_eventN1_tIIw_bdt", "Highest Disc Cand (Type IIw) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_IImiq_bdt = new TH2F ("h_eventN1_tIImiq_bdt", "Highest Disc Cand (Type IImiq) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_Itp_bdt = new TH2F ("h_eventN1_tItp_bdt", "Highest Disc Cand (Type It / Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_It_bdt = new TH2F ("h_eventN1_tIt_bdt", "Highest Disc Cand (Type It) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_Ip_bdt = new TH2F ("h_eventN1_tIp_bdt", "Highest Disc Cand (Type Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_0xI_bdt = new TH2F ("h_eventN1_t0x1_bdt", "Highest Disc Cand (Type 0 / I) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_I_bdt = new TH2F ("h_eventN1_t1_bdt", "Highest Disc Cand (Type 1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN1_0_bdt = new TH2F ("h_eventN1_t0_bdt", "Highest Disc Cand (Type 0) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);


  TH2F *h_eventN2_sum_bdt = new TH2F ("h_eventN2_sum_bdt", "2nd Highest Disc Cand (all) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_III_bdt = new TH2F ("h_eventN2_tIII_bdt", "2nd Highest Disc Cand (Type III) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_II_bdt = new TH2F ("h_eventN2_tII_bdt", "2nd Highest Disc Cand (Type II) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_IIb_bdt = new TH2F ("h_eventN2_tIIb_bdt", "2nd Highest Disc Cand (Type IIb) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_IImib_bdt = new TH2F ("h_eventN2_tIImib_bdt", "2nd Highest Disc Cand (Type IImib) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_IIw_bdt = new TH2F ("h_eventN2_tIIw_bdt", "2nd Highest Disc Cand (Type IIw) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_IImiq_bdt = new TH2F ("h_eventN2_tIImiq_bdt", "2nd Highest Disc Cand (Type IImiq) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_Itp_bdt = new TH2F ("h_eventN2_tItp_bdt", "2nd Highest Disc Cand (Type It / Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_It_bdt = new TH2F ("h_eventN2_tIt_bdt", "2nd Highest Disc Cand (Type It) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_Ip_bdt = new TH2F ("h_eventN2_tIp_bdt", "2nd Highest Disc Cand (Type Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_0xI_bdt = new TH2F ("h_eventN2_t0x1_bdt", "2nd Highest Disc Cand (Type 0 / 1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_I_bdt = new TH2F ("h_eventN2_t1_bdt", "2nd Highest Disc Cand (Type 1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN2_0_bdt = new TH2F ("h_eventN2_t0_bdt", "2nd Highest Disc Cand (Type 0) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);

  TH2F *h_eventN3_sum_bdt = new TH2F ("h_eventN3_sum_bdt", "3rd Highest Disc Cand (all) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_III_bdt = new TH2F ("h_eventN3_tIII_bdt", "3rd Highest Disc Cand (Type III) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_II_bdt = new TH2F ("h_eventN3_tII_bdt", "3rd Highest Disc Cand (Type II) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_IIb_bdt = new TH2F ("h_eventN3_tIIb_bdt", "3rd Highest Disc Cand (Type IIb) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_IImib_bdt = new TH2F ("h_eventN3_tIImib_bdt", "3rd Highest Disc Cand (Type IImib) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_IIw_bdt = new TH2F ("h_eventN3_tIIw_bdt", "3rd Highest Disc Cand (Type IIw) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_IImiq_bdt = new TH2F ("h_eventN3_tIImiq_bdt", "3rd Highest Disc Cand (Type IImiq) ;Discriminant; Number of Jets", 
				     nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_Itp_bdt = new TH2F ("h_eventN3_tItp_bdt", "3rd Highest Disc Cand (Type It / Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_It_bdt = new TH2F ("h_eventN3_tIt_bdt", "3rd Highest Disc Cand (Type It) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_Ip_bdt = new TH2F ("h_eventN3_tIp_bdt", "3rd Highest Disc Cand (Type Ip) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_0xI_bdt = new TH2F ("h_eventN3_t0x1_bdt", "3rd Highest Disc Cand (Type 0x1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_I_bdt = new TH2F ("h_eventN3_t1_bdt", "3rd Highest Disc Cand (Type 1) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);
  TH2F *h_eventN3_0_bdt = new TH2F ("h_eventN3_t0_bdt", "3rd Highest Disc Cand (Type 0) ;Discriminant; Number of Jets", 
				      nBins, 0.0, 1.0, 7, 7, 14);

  //Count number of correct hadronic top quark jets that are tagged, by jet multiplicity and discriminant
  TH3F *h_eventN1_CntJet_hot = new TH3F ("h_eventN1_CntJet_hot", "N1 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 4, 0, 4);
  TH3F *h_eventN2_CntJet_hot = new TH3F ("h_eventN2_CntJet_hot", "N2 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 4, 0, 4);
  TH3F *h_eventN3_CntJet_hot = new TH3F ("h_eventN3_CntJet_hot", "N3 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 4, 0, 4);
  TH3F *h_eventN12_CntJet_hot = new TH3F ("h_eventN12_CntJet_hot", "N1 + N2 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 7, 0, 7);
  TH3F *h_eventN123_CntJet_hot = new TH3F ("h_eventN123_CntJet_hot", "N1 + N2 + N3 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 10, 0, 10);

  TH3F *h_eventN1_CntJet_bdt = new TH3F ("h_eventN1_CntJet_bdt", "N1 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 4, 0, 4);
  TH3F *h_eventN2_CntJet_bdt = new TH3F ("h_eventN2_CntJet_bdt", "N2 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 4, 0, 4);
  TH3F *h_eventN3_CntJet_bdt = new TH3F ("h_eventN3_CntJet_bdt", "N3 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 4, 0, 4);
  TH3F *h_eventN12_CntJet_bdt = new TH3F ("h_eventN12_CntJet_bdt", "N1 + N2 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 7, 0, 7);
  TH3F *h_eventN123_CntJet_bdt = new TH3F ("h_eventN123_CntJet_bdt", "N1 + N2 + N3 Aggregate;Discriminant; Number of Jets in Event; Aggregate of Tagged Hadronic Top Jets", 
					 nBins, 0.0, 1.0, 7, 7, 14, 10, 0, 10);
  
  //BDT Scaling histograms (map range of BDT to [0,1] using function scalebdt(BDTdiscriminant)
  TH2F *h_noscaling_bdt = new TH2F ("h_noscaling_bdt", "Un-scaled BDT Distribution;Discriminant;Number of Tagger Candidates", 200, -0.45, 0.35, 7, 7, 14);
  TH2F *h_noscaling2_bdt = new TH2F ("h_noscaling2_bdt", "Un-scaled BDT Distribution (2nd Highest);Discriminant;Number of Tagger Candidates", 200, -0.45, 0.35, 7, 7, 14);
  TH2F *h_scaling_bdt = new TH2F ("h_scaling_bdt", "scaled BDT Distribution;Discriminant;Number of Tagger Candidates", 200, -0.2, 1.4, 7, 7, 14);
  TH2F *h_scaling2_bdt = new TH2F ("h_scaling2_bdt", "scaled BDT Distribution (2nd Highest);Discriminant;Number of Tagger Candidates", 200, 0.0, 1.0, 7, 7, 14);

  //store errors
  TH1I *h_nERROR = new TH1I ("h_nERROR", "Number of logic errors in classification system + miscellaneous", 100, 0, 100);

  //TH1F *TypeWDiscr = new TH1F ("h_fullWdiscr", "Fully Reconstructible W bosons; Discriminant; Number of matches", nBins, 0.0, 1.0);

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
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/200kTwoTop.root");

    //New style files
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/NewTT50HT3J.root");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/NewTT50HT.root");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/NewTT250HT.root");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/NewTT500HT.root");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/NewTTTT500HT.root");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/ShortTT500HT_corr.root", "r"); //test of corrections to b-tagging, eta, pt requirements
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/StoreTT50HT.root", "r");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/StoreTT500HT.root", "r");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/StoreTTSLFilt500HT.root", "r");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/StoreTTTT500HT.root", "r");

    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/LongTT500HT_corr.root", "r"); //corrected btagging sample
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/LongTTSLFilt500HT_corr.root", "r");
    //tf2 = TFile::Open("/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/LongTTTT500HT_corr.root", "r");

    //Get TFile from input string
    tf2 = TFile::Open(inFile.c_str(), "r");

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

    //Custom variables for reading newer style gen-matched jets and status flags
    std::vector<TLorentzVector>** recoTop1 = new std::vector<TLorentzVector>*();
    std::vector<TLorentzVector>** recoTop2 = new std::vector<TLorentzVector>*();
    std::vector<TLorentzVector>** recoTop3 = new std::vector<TLorentzVector>*();
    std::vector<TLorentzVector>** recoTop4 = new std::vector<TLorentzVector>*();
    std::vector<int>** recoTop1flags = new std::vector<int>*();
    std::vector<int>** recoTop2flags = new std::vector<int>*();
    std::vector<int>** recoTop3flags = new std::vector<int>*();
    std::vector<int>** recoTop4flags = new std::vector<int>*();

    //run, lumi, event variables
    uint run;
    uint lumiBlock;
    uint event;

    //Deactivate all branches, then activate the branches of interest
    tree->SetBranchStatus("*", 0);
    std::cout << "Deactivated all branches" << std::endl;

    //standard exampleInputs.root settings
    if(runStdExample == true){
      //Activate branches of interest
      //AK4 jet lorentz vectors
      tree->SetBranchStatus( "ak4jetsLVec", 1);
      tree->SetBranchAddress("ak4jetsLVec", AK4JetLV);
    
      //AK4 jet b-tag values
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
      //run, lumi, event branches
      tree->SetBranchStatus( "run", 1);
      tree->SetBranchAddress("run", &run);

      tree->SetBranchStatus( "lumiBlock", 1);
      tree->SetBranchAddress("lumiBlock", &lumiBlock);

      tree->SetBranchStatus( "event", 1);
      tree->SetBranchAddress("event", &event);

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

      if(isNewStyle){
	//New style tops' TLV's
	//top 1
	tree->SetBranchStatus( "recoTop1", 1);
	tree->SetBranchAddress("recoTop1", recoTop1);
	//top 2
	tree->SetBranchStatus( "recoTop2", 1);
	tree->SetBranchAddress("recoTop2", recoTop2);
	//top 3
	tree->SetBranchStatus( "recoTop3", 1);
	tree->SetBranchAddress("recoTop3", recoTop3);
	//top 4
	tree->SetBranchStatus( "recoTop4", 1);
	tree->SetBranchAddress("recoTop4", recoTop4);

	//top 1
	tree->SetBranchStatus( "recoTop1flags", 1);
	tree->SetBranchAddress("recoTop1flags", recoTop1flags);
	//top 2
	tree->SetBranchStatus( "recoTop2flags", 1);
	tree->SetBranchAddress("recoTop2flags", recoTop2flags);
	//top 3
	tree->SetBranchStatus( "recoTop3flags", 1);
	tree->SetBranchAddress("recoTop3flags", recoTop3flags);
	//top 4
	tree->SetBranchStatus( "recoTop4flags", 1);
	tree->SetBranchAddress("recoTop4flags", recoTop4flags);
      }

      std::cout << "tree branches attached" << std::endl;
    }

    //Create top tagger object
    TopTagger tt;

    //Setup TMVA for evaluating the 2016 default BDT
    float defaultVars[6] = {0};
    TMVA::Reader *reader;
    reader = new TMVA::Reader( "!Color:!Silent" );
    reader->AddVariable("btag", &defaultVars[0]);
    reader->AddVariable("ThPtOverSumPt", &defaultVars[1]);
    reader->AddVariable("AngleThWh", &defaultVars[2]);
    reader->AddVariable("AngleThBh", &defaultVars[3]);
    reader->AddVariable("HadrWmass", &defaultVars[4]);
    reader->AddVariable("TopMass", &defaultVars[5]);
    reader->BookMVA( "BDT method", "JetCombTrainer_BDT.weights.xml" );
    //Reader is ready to evaluate.
    defaultVars[0] = 0.99;
    defaultVars[0] = 0.6;
    defaultVars[0] = 0.13;
    defaultVars[0] = 0.3;
    defaultVars[0] = 81.0;
    defaultVars[0] = 173.2;

    std::cout << "Testing MVA: " << reader->EvaluateMVA( "BDT method" ) << std::endl;
    
    //try-catch on TTException which are thrown by the top tagger
    try
    {
        //Set top tagger cfg file
        tt.setCfgFile("TopTagger.cfg");

	//debug counter
	int nERROR = 0;

        //Loop over events
        int Nevt = 0;
        while(tree->GetEntry(Nevt))
        {
            //increment event number
            ++Nevt;
	    //Fixed buggy logic comparison with -1, switching to != -1 statement
	    if(maxEventsToProcess != -1 && maxEventsToProcess < Nevt)
	      break;

            //Print event number 
            if(!silentRunning)
	      printf("\n===================================================>\nEvent #: %i\n", Nevt);

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
		
	    //count reconstructible tops via the flags and other means...
	    int nAnyHadTops = 0;
	    int nRecoTops = 0;
	    int nHOTTops = 0;
	    int nBDTTops = 0;

	    //get number of jets for categorization of discriminants
	    int nJets = (**AK4JetLV).size();
	    h_nJet->Fill(nJets);


	    //Create ResTTPermuter, be aware of the Most Vexing Parse (C++), necessitating curly brackets around argument if it looks like a function declaration
	    //Create loopable version
	    ResTTPermuter Pitcher(*AK4JetLV, *AK4JetBtag);
	    uint loopLimit = Pitcher.GetNumTrijets();
	    for(uint loop = 0; loop < loopLimit; loop++){
	      //std::cout << ">>>>>>>>>>>>>>>>>>>>> Loop: " << loop << std::endl;
	      std::vector<BDTCand> lineup = Pitcher.PitchCandidates(loop);
	      for(int bat = 0; bat < lineup.size(); bat++){
		if(debug2){
		  std::cout << "\t>>"
			    << lineup[bat].idx_i << " "
			    << lineup[bat].idx_j << " "
			    << lineup[bat].idx_k << " "
			    << lineup[bat].btag << " "
			    << lineup[bat].ThPtOverSumPt << " "
			    << lineup[bat].AngleThWh << " "
			    << lineup[bat].AngleThBh << " "
			    << lineup[bat].HadrWmass << " "
			    << lineup[bat].TopMass << " "
			    << lineup[bat].discriminant << " ";
		}
		//<< std::endl;
		defaultVars[0] = lineup[bat].btag;
		defaultVars[1] = lineup[bat].ThPtOverSumPt;
		defaultVars[2] = lineup[bat].AngleThWh;
		defaultVars[3] = lineup[bat].AngleThBh;
		defaultVars[4] = lineup[bat].HadrWmass;
		defaultVars[5] = lineup[bat].TopMass;
		lineup[bat].discriminant = scaleBDT(reader->EvaluateMVA( "BDT method" ));
		if(debug2)
		  std::cout << "   -->   " << lineup[bat].discriminant << std::endl;
	      }
	      std::sort(lineup.begin(), lineup.end(), [](const BDTCand &left, const BDTCand &right) {
		  return left.discriminant > right.discriminant;
		});
	      // for(int bat = 0; bat < lineup.size(); bat++){
	      // 	std::cout << "lineup disc: " << lineup[bat].discriminant << std::endl;
	      // }
	      Pitcher.CatchCandidate(lineup[0]);
	    }
	    auto BDTResults = Pitcher.GetFinalCandidates();
	    nBDTTops = BDTResults.size();
	    if(nBDTTops != Pitcher.GetNumTrijets())
	      nERROR++;
	    
	    

	    //no longer really used...
	    int nRecoWs = 0;
	    bool Top1Reco = false;
	    bool Top2Reco = false;
	    bool Top3Reco = false;
	    if(debug1) printf("Booleans for Top reconstruction (initialization) %2d %2d %2d", Top1Reco, Top2Reco, Top3Reco);
	    uint flagcounter = 0;
	    if(debug1) printf("\n\tDebugging reconstructible top counting and flags: ");


	    // //Testing the matching types
	    // TLorentzVector b1, b2, b3, q11, q12, q13, q21, q22, q23, oth;
	    // b1.SetPtEtaPhiE(45, 1.234, -2.14, 48);
	    // b2.SetPtEtaPhiE(56, 2.341, 2.44, 58.5);
	    // b3.SetPtEtaPhiE(32, -1.154, -0.19, 37);
	    // q11.SetPtEtaPhiE(33, 0.154, 0.719, 35);
	    // q12.SetPtEtaPhiE(90, 0.548, -0.103, 99);
	    // q13.SetPtEtaPhiE(63, -2.17, 0.356, 70);
	    // q21.SetPtEtaPhiE(23, -0.06, 0.056, 25);
	    // q22.SetPtEtaPhiE(80, 0.958, 0.561, 80.99);
	    // q23.SetPtEtaPhiE(105, -0.481, -1.832, 110);
	    // oth.SetPtEtaPhiE(120, 0.0, 1.23, 121);
	    // std::vector<TLorentzVector> *Reco1, *Reco2, *Reco3;
	    // std::vector<TLorentzVector> *tIII, *tIIb, *tIIw;
	    // std::vector<TLorentzVector> *tIImib, *tIImiq, *tI, *t0;
	    // std::vector<int> bling;
	    // bling.push_back(2);
	    // Reco1 = new std::vector<TLorentzVector>();
	    // Reco2 = new std::vector<TLorentzVector>();
	    // Reco3 = new std::vector<TLorentzVector>();
	    // tIII = new std::vector<TLorentzVector>();
	    // tIIb = new std::vector<TLorentzVector>();
	    // tIIw = new std::vector<TLorentzVector>();
	    // tIImib = new std::vector<TLorentzVector>();
	    // tIImiq = new std::vector<TLorentzVector>();
	    // tI = new std::vector<TLorentzVector>();
	    // t0 = new std::vector<TLorentzVector>();
	    // Reco1->push_back(b1);
	    // Reco1->push_back(q11);
	    // Reco1->push_back(q21);

	    // Reco2->push_back(b2);
	    // Reco2->push_back(q12);
	    // Reco2->push_back(q22);

	    // Reco3->push_back(b3);
	    // Reco3->push_back(q13);
	    // Reco3->push_back(q23);

	    // tIII->push_back(b1);
	    // tIII->push_back(q11);
	    // tIII->push_back(q21);

	    // tIIb->push_back(b3);
	    // tIIb->push_back(q12);
	    // tIIb->push_back(q22);
	
	    // tIImib->push_back(oth);
	    // tIImib->push_back(q13);
	    // tIImib->push_back(q23);	

	    // tIIw->push_back(b2);
	    // tIIw->push_back(q11);
	    // tIIw->push_back(q22);
	
	    // tIImiq->push_back(b3);
	    // tIImiq->push_back(q23);
	    // tIImiq->push_back(oth);	

	    // tI->push_back(b1);
	    // tI->push_back(q12);
	    // tI->push_back(q23);
	
	    // t0->push_back(oth);
	    // t0->push_back(oth);
	    // t0->push_back(oth);	

	    // ResTTEvaluator Testing("Testing", true, false);
	    // Testing.addReco(Reco1, bling);
	    // Testing.addReco(Reco2, bling);
	    // Testing.addReco(Reco3, bling);
	    // Testing.addCand(tIII, 0.95);
	    // Testing.addCand(tIIb, 0.90);
	    // Testing.addCand(tIIw, 0.85);
	    // Testing.addCand(tIImib, 0.80);
	    // Testing.addCand(tIImiq, 0.75);
	    // Testing.addCand(tI, 0.70);
	    // Testing.addCand(t0, 0.65);

	    // Testing.evaluateR();


	    //ResTTEvaluator HOTEval("HOT");
	    ResTTEvaluator HOTEval("HOT", false, false);
	    
	    ResTTEvaluator BDTEval("BDT", false, false);

	    //Old style tops
	    if(!isNewStyle){
	      std::vector<int> Top1flags;
	      Top1flags.push_back(0);
	      auto the1Top = **hadTop1Constit;
	      if(the1Top.size() > 0){
		HOTEval.addReco(*hadTop1Constit, Top1flags);
		BDTEval.addReco(*hadTop1Constit, Top1flags);
		//nAnyHadTops++;
	      }
	      auto the2Top = **hadTop2Constit;
	      if(the2Top.size() > 0){
		HOTEval.addReco(*hadTop2Constit, Top1flags);
		BDTEval.addReco(*hadTop2Constit, Top1flags);
		//nAnyHadTops++;
	      }
	      auto the3Top = **hadTop3Constit;
	      if(the3Top.size() > 0){
		HOTEval.addReco(*hadTop3Constit, Top1flags);
		BDTEval.addReco(*hadTop3Constit, Top1flags);
		//nAnyHadTops++;
	      }
	    }

	    //new style tops
	    if(isNewStyle){
	      if((**recoTop1).size() > 0 && (**recoTop1flags)[0] > -1){
		//if((**recoTop1).size() > 0){ leptonic permitted
		HOTEval.addReco(*recoTop1, **recoTop1flags);
		BDTEval.addReco(*recoTop1, **recoTop1flags);
		if((**recoTop1flags)[0] > -1){
		  nAnyHadTops++;
		  if((**recoTop1flags)[0] > 2)
		    nRecoTops++;
		}
	      }
	      if((**recoTop2).size() > 0 && (**recoTop2flags)[0] > -1){
		//if((**recoTop2).size() > 0){
		HOTEval.addReco(*recoTop2, **recoTop2flags);
		BDTEval.addReco(*recoTop2, **recoTop2flags);
		if((**recoTop2flags)[0] > -1){
		  nAnyHadTops++;
		  if((**recoTop2flags)[0] > 2)
		    nRecoTops++;
		}
	      }
	      if((**recoTop3).size() > 0 && (**recoTop3flags)[0] > -1){
		//if((**recoTop3).size() > 0){
		HOTEval.addReco(*recoTop3, **recoTop3flags);
		BDTEval.addReco(*recoTop3, **recoTop3flags);
		if((**recoTop3flags)[0] > -1){
		  nAnyHadTops++;
		  if((**recoTop3flags)[0] > 2)
		    nRecoTops++;
		}
	      }
	      if((**recoTop4).size() > 0 && (**recoTop4flags)[0] > -1){
		//if((**recoTop4).size() > 0){
		HOTEval.addReco(*recoTop4, **recoTop4flags);
		BDTEval.addReco(*recoTop4, **recoTop4flags);
		if((**recoTop4flags)[0] > -1){
		  nAnyHadTops++;
		  if((**recoTop4flags)[0] > 2)
		    nRecoTops++;
		}
	      }
	      if(nAnyHadTops == 0){
		if(debug1)
		  std::cout << "\n>>>>>>>>>>>>\n";
		for(const uint oldflag: **FlagTop)
		  std::cout << "\tO: " << oldflag;
		if(debug1)
		  std::cout << "\n\tRun: " << run << "\tLumi: " << lumiBlock << "\tEvent: " << event << std::endl;
	      }
	    }
	    if(!silentRunning)
	      std::cout << "\n\t\t\t\tNumber of Hadronic Tops passed from event: " << nAnyHadTops << std::endl;

	    // if(the1Top.size() > 0)
	    //   HOTEval.addCand(*hadTop1Constit, 0.873);
	    // if(the2Top.size() > 0)
	    //   HOTEval.addCand(*hadTop2Constit, 0.233);
	    // if(the3Top.size() > 0)
	    //   HOTEval.addCand(*hadTop3Constit, 0.463);

	    if(!silentRunning)
	      std::cout << "\n============================" << std::endl;
	    for(const uint topflag: **FlagTop){
	      if(debug1) printf("\n\ttopflag = %6d", topflag);
	      if(topflag < 9999){
		if(topflag > 1020){
		  RecoTypes->Fill(2);
		  //nRecoTops += 1;
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
	    //printf("\n\t\t\t# Fully Reconstuctible Tops: %2d \t Reconstructible W: %2d\n", nRecoTops, nRecoWs);

	    
	    nHOTTops = 0;
            //print top properties
            for(const TopObject* top : tops)
            {
	        nHOTTops++;
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

		//prepare container for candidate:
		if(debug2)
		  std::cout << "=========> Top candidate" << std::endl;
		std::vector<TLorentzVector>* theCand;
		theCand = new std::vector<TLorentzVector>;

                //Print properties of individual top constituent jets 
                for(const Constituent* constituent : constituents)
                {
		  //package jets in the candidate
		  theCand->push_back(constituent->p());
		  
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
		      mat2 += (inJet == constituent->p());
		    }
		    //printf("\n");
		    for(const TLorentzVector inJet : **hadTop2Constit){
		      mat3 += (inJet == constituent->p());
		    }
		    //printf("\n");
		    for(const TLorentzVector inJet : **hadTop3Constit){
		      mat4 += (inJet == constituent->p());
		    }
                }

		//add candidate to ResTTEvaluator for HOT
		HOTEval.addCand(theCand, top->getDiscriminator());
		//printf("\t\t\ttop1: %2d || top2: %2d|| top3: %2d || Discriminant: %6.4lf\n", mat2, mat3, mat4, top->getDiscriminator());
		//if(mat2 == 3 || mat3 == 3 || mat4 == 3) h_typeIII->Fill(top->getDiscriminator());
		//if(mat2 == 2 || mat3 == 2 || mat4 == 2) TypeWDiscr->Fill(top->getDiscriminator());
		//if(mat2 < 2 && mat3 < 2 && mat4 < 2) h_typeI->Fill(top->getDiscriminator());
		delete theCand;
            }
	    
	    //add candidate to ResTTEvaluator for BDT, results were stored in BDTResults as vector of pairs of TLVs and discriminant
	    for(uint bdtrescnt = 0; bdtrescnt < BDTResults.size(); bdtrescnt++)
	      BDTEval.addCand(&BDTResults[bdtrescnt].first, BDTResults[bdtrescnt].second);
	    if(debug2)
	      BDTEval.printDimensions();
	    BDTEval.evaluateR();
	    std::vector<std::pair<double, std::string>> tCR_bdt = BDTEval.getOrderedClassesR();
	    std::vector<std::pair<double, int>> tAR_bdt = BDTEval.getOrderedAggregatesR();
	    //*

	    if(tAR_bdt.size() > 0){
	      h_eventN1_CntJet_bdt->Fill(tAR_bdt[0].first, nJets, tAR_bdt[0].second);
	    }
	    if(tAR_bdt.size() > 1){
	      h_eventN2_CntJet_bdt->Fill(tAR_bdt[1].first, nJets, tAR_bdt[1].second);
	      h_eventN12_CntJet_bdt->Fill(tAR_bdt[1].first, nJets, tAR_bdt[0].second + tAR_bdt[1].second);
	    }
	    if(tAR_bdt.size() > 2){
	      h_eventN3_CntJet_bdt->Fill(tAR_bdt[2].first, nJets, tAR_bdt[2].second);
	      h_eventN123_CntJet_bdt->Fill(tAR_bdt[2].first, nJets, tAR_bdt[0].second + tAR_bdt[1].second + tAR_bdt[2].second);
	    }

	    for(int e = 0; e < tCR_bdt.size(); e++){
	      h_scaling_bdt->Fill(tCR_bdt[e].first, nJets);
	      h_noscaling_bdt->Fill( unscaleBDT(tCR_bdt[e].first), nJets );
	      if(tCR_bdt[e].second == "typeIII")
		h_typeIII_bdt->Fill(tCR_bdt[e].first, nJets);
	      else if(tCR_bdt[e].second == "typeIIb")
		h_typeIIb_bdt->Fill(tCR_bdt[e].first, nJets);
	      else if(tCR_bdt[e].second == "typeIImib")
		h_typeIImib_bdt->Fill(tCR_bdt[e].first, nJets);
	      else if(tCR_bdt[e].second == "typeIIw")
		h_typeIIw_bdt->Fill(tCR_bdt[e].first, nJets);
	      else if(tCR_bdt[e].second == "typeIImiq")
		h_typeIImiq_bdt->Fill(tCR_bdt[e].first, nJets);
	      else if(tCR_bdt[e].second == "typeIt")
		h_typeIt_bdt->Fill(tCR_bdt[e].first, nJets);
	      else if(tCR_bdt[e].second == "typeIp")
		h_typeIp_bdt->Fill(tCR_bdt[e].first, nJets);
	      else if(tCR_bdt[e].second == "typeI" || tCR_bdt[e].second == "type0")
		h_type0xI_bdt->Fill(tCR_bdt[e].first, nJets);
	      else{
		std::cout << "LOGIC FAILURE! LOGIC FAILURE! LOGIC FAILURE!" << std::endl;
		std::cout << tCR_bdt[e].second << std::endl;
		nERROR++;
		  }
	      //ordered candidates
	      switch(e){
	      case 0:{
		h_nJet_N1_bdt->Fill(nJets);
		if(tCR_bdt[e].second == "typeIII"){
		  h_eventN1_III_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIIb"){
		  h_eventN1_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_IIb_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIImib"){
		  h_eventN1_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_IImib_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIIw"){
		  h_eventN1_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_IIw_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIImiq"){
		  h_eventN1_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_IImiq_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIt"){
		  h_eventN1_Itp_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_It_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIp"){
		  h_eventN1_Itp_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_Ip_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeI"){
		  h_eventN1_0xI_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_I_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "type0"){
		  h_eventN1_0xI_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_0_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN1_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else{
		  std::cout << "LOGIC FAILURE! LOGIC FAILURE! LOGIC FAILURE!" << std::endl;
		  std::cout << tCR_bdt[e].second << std::endl;
		  nERROR++;
		}
	      }//case 0 (Highest)
		break;
	      case 1:{
		h_nJet_N2_bdt->Fill(nJets);
		h_scaling2_bdt->Fill(tCR_bdt[e].first, nJets);
		h_noscaling2_bdt->Fill( unscaleBDT(tCR_bdt[e].first), nJets );
		if(tCR_bdt[e].second == "typeIII"){
		  h_eventN2_III_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIIb"){
		  h_eventN2_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_IIb_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIImib"){
		  h_eventN2_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_IImib_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIIw"){
		  h_eventN2_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_IIw_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIImiq"){
		  h_eventN2_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_IImiq_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIt"){
		  h_eventN2_Itp_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_It_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIp"){
		  h_eventN2_Itp_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_Ip_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeI"){
		  h_eventN2_0xI_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_I_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "type0"){
		  h_eventN2_0xI_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_0_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN2_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else{
		  std::cout << "LOGIC FAILURE! LOGIC FAILURE! LOGIC FAILURE!" << std::endl;
		  std::cout << tCR_bdt[e].second << std::endl;
		  nERROR++;
		}
	      }//case 1 (2nd highest)
		break;
	      case 2:{
		h_nJet_N3_bdt->Fill(nJets);
		if(tCR_bdt[e].second == "typeIII"){
		  h_eventN3_III_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIIb"){
		  h_eventN3_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_IIb_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIImib"){
		  h_eventN3_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_IImib_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIIw"){
		  h_eventN3_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_IIw_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIImiq"){
		  h_eventN3_II_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_IImiq_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIt"){
		  h_eventN3_Itp_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_It_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeIp"){
		  h_eventN3_Itp_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_Ip_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "typeI"){
		  h_eventN3_0xI_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_I_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else if(tCR_bdt[e].second == "type0"){
		  h_eventN3_0xI_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_0_bdt->Fill(tCR_bdt[e].first, nJets);
		  h_eventN3_sum_bdt->Fill(tCR_bdt[e].first, nJets);
		}
		else{
		  std::cout << "LOGIC FAILURE! LOGIC FAILURE! LOGIC FAILURE!" << std::endl;
		  std::cout << tCR_bdt[e].second << std::endl;
		  nERROR++;
		}
	      }//end case 2 (3rd highest)
		break;
	      default: break;
	      }
	    } //*/


	    if(debug2){
	      std::cout << "\nNumber of HOT Tagger Candidates: " << nHOTTops << std::endl;
	      HOTEval.printDimensions();
	    }
	    HOTEval.evaluateR();
	    std::vector<std::pair<double, std::string>> tCR_hot = HOTEval.getOrderedClassesR();
	    std::vector<std::pair<double, int>> tAR_hot = HOTEval.getOrderedAggregatesR();
	    
	    //Error checking
	    if(debug2)
	      std::cout << "\n\tOrdered candidates: ";
	    if(tCR_hot.size() > 1)
	      for(int f = 2; f < tCR_hot.size(); f++){
		//std::cout << tCR_hot[f].first << " (" << tCR_hot[f].second << ")\t";
		if(tCR_hot[f].first > tCR_hot[f-1].first) nERROR++;
	      }
	    if(tCR_hot.size() != nHOTTops) nERROR++;

	    //Fill Aggregate Jet histograms according to discriminant, nJets, and aggregate number of matched hadronic top jets at that discr. point and above
	    if(tAR_hot.size() > 0){
	      h_eventN1_CntJet_hot->Fill(tAR_hot[0].first, nJets, tAR_hot[0].second);
	    }
	    if(tAR_hot.size() > 1){
	      h_eventN2_CntJet_hot->Fill(tAR_hot[1].first, nJets, tAR_hot[1].second);
	      h_eventN12_CntJet_hot->Fill(tAR_hot[1].first, nJets, tAR_hot[0].second + tAR_hot[1].second);
	    }
	    if(tAR_hot.size() > 2){
	      h_eventN3_CntJet_hot->Fill(tAR_hot[2].first, nJets, tAR_hot[2].second);
	      h_eventN123_CntJet_hot->Fill(tAR_hot[2].first, nJets, tAR_hot[0].second + tAR_hot[1].second + tAR_hot[2].second);
	    }

	    for(int e = 0; e < tCR_hot.size(); e++){
	      if(tCR_hot[e].second == "typeIII")
		h_typeIII_hot->Fill(tCR_hot[e].first, nJets);
	      else if(tCR_hot[e].second == "typeIIb")
		h_typeIIb_hot->Fill(tCR_hot[e].first, nJets);
	      else if(tCR_hot[e].second == "typeIImib")
		h_typeIImib_hot->Fill(tCR_hot[e].first, nJets);
	      else if(tCR_hot[e].second == "typeIIw")
		h_typeIIw_hot->Fill(tCR_hot[e].first, nJets);
	      else if(tCR_hot[e].second == "typeIImiq")
		h_typeIImiq_hot->Fill(tCR_hot[e].first, nJets);
	      else if(tCR_hot[e].second == "typeIt")
		h_typeIt_hot->Fill(tCR_hot[e].first, nJets);
	      else if(tCR_hot[e].second == "typeIp")
		h_typeIp_hot->Fill(tCR_hot[e].first, nJets);
	      else if(tCR_hot[e].second == "typeI" || tCR_hot[e].second == "type0")
		h_type0xI_hot->Fill(tCR_hot[e].first, nJets);
	      else{
		std::cout << "LOGIC FAILURE! LOGIC FAILURE! LOGIC FAILURE!" << std::endl;
		std::cout << tCR_hot[e].second << std::endl;
		nERROR++;
		  }
	      //ordered candidates
	      switch(e){
	      case 0:{
		h_nJet_N1_hot->Fill(nJets);
		if(tCR_hot[e].second == "typeIII"){
		  h_eventN1_III_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIIb"){
		  h_eventN1_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_IIb_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIImib"){
		  h_eventN1_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_IImib_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIIw"){
		  h_eventN1_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_IIw_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIImiq"){
		  h_eventN1_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_IImiq_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIt"){
		  h_eventN1_Itp_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_It_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIp"){
		  h_eventN1_Itp_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_Ip_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeI"){
		  h_eventN1_0xI_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_I_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "type0"){
		  h_eventN1_0xI_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_0_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN1_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else{
		  std::cout << "LOGIC FAILURE! LOGIC FAILURE! LOGIC FAILURE!" << std::endl;
		  std::cout << tCR_hot[e].second << std::endl;
		  nERROR++;
		}
	      }//case 0 (Highest)
		break;
	      case 1:{
		h_nJet_N2_hot->Fill(nJets);
		if(tCR_hot[e].second == "typeIII"){
		  h_eventN2_III_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIIb"){
		  h_eventN2_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_IIb_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIImib"){
		  h_eventN2_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_IImib_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIIw"){
		  h_eventN2_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_IIw_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIImiq"){
		  h_eventN2_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_IImiq_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIt"){
		  h_eventN2_Itp_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_It_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIp"){
		  h_eventN2_Itp_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_Ip_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeI"){
		  h_eventN2_0xI_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_I_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "type0"){
		  h_eventN2_0xI_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_0_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN2_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else{
		  std::cout << "LOGIC FAILURE! LOGIC FAILURE! LOGIC FAILURE!" << std::endl;
		  std::cout << tCR_hot[e].second << std::endl;
		  nERROR++;
		}
	      }//case 1 (2nd highest)
		break;
	      case 2:{
		h_nJet_N3_hot->Fill(nJets);
		if(tCR_hot[e].second == "typeIII"){
		  h_eventN3_III_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIIb"){
		  h_eventN3_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_IIb_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIImib"){
		  h_eventN3_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_IImib_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIIw"){
		  h_eventN3_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_IIw_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIImiq"){
		  h_eventN3_II_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_IImiq_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIt"){
		  h_eventN3_Itp_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_It_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeIp"){
		  h_eventN3_Itp_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_Ip_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "typeI"){
		  h_eventN3_0xI_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_I_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else if(tCR_hot[e].second == "type0"){
		  h_eventN3_0xI_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_0_hot->Fill(tCR_hot[e].first, nJets);
		  h_eventN3_sum_hot->Fill(tCR_hot[e].first, nJets);
		}
		else{
		  std::cout << "LOGIC FAILURE! LOGIC FAILURE! LOGIC FAILURE!" << std::endl;
		  std::cout << tCR_hot[e].second << std::endl;
		  nERROR++;
		}
	      }//end case 2 (3rd highest)
		break;
	      default: break;
	      }
	    }

            //Print properties of the remaining system
            //the remaining system is used as the second portion of the visible system to calculate MT2 in the NT = 1 bin
            //const TopObject& rsys = ttr.getRsys();
            //printf("\tRsys properties: N constituents: %3d,   Pt: %6.1lf,   Eta: %7.3lf,   Phi: %7.3lf\n", rsys.getNConstituents(), rsys.p().Pt(), rsys.p().Eta(), rsys.p().Phi());
        
	    //Fill Event Histos here
	    h_nTrueRecoTops_full->Fill(nRecoTops);
	    h_nTrueRecoTops->Fill(nAnyHadTops);
	    h_nCandTops_hot->Fill(nHOTTops);
	    h_nCandTops_bdt->Fill(nBDTTops);
	    switch (nAnyHadTops) {
	    case 1: 
	      h_nCandTops_hot1->Fill(nHOTTops);
	      h_nCandTops_bdt1->Fill(nBDTTops);
	      break;
	    case 2: 
	      h_nCandTops_hot2->Fill(nHOTTops);
	      h_nCandTops_bdt2->Fill(nBDTTops);
	      break;
	    case 3: 
	      h_nCandTops_hot3->Fill(nHOTTops);
	      h_nCandTops_bdt3->Fill(nBDTTops);
	      break;
	    case 4: 
	      h_nCandTops_hot4->Fill(nHOTTops);
	      h_nCandTops_bdt4->Fill(nBDTTops);
	      break;
	    default:
	      h_nCandTops_hot_oth->Fill(nHOTTops);
	      h_nCandTops_bdt_oth->Fill(nBDTTops);
	      break;
	    }
	    switch (nRecoTops) { //only fully reconstructable tops... not good for ensuring everything adds up!
	    case 0:
	      h_nCandFullTops_hot_0->Fill(nHOTTops);
	      h_nCandFullTops_bdt_0->Fill(nBDTTops);
	      break;
	    case 1: 
	      h_nCandFullTops_hot1->Fill(nHOTTops);
	      h_nCandFullTops_bdt1->Fill(nBDTTops);
	      break;
	    case 2: 
	      h_nCandFullTops_hot2->Fill(nHOTTops);
	      h_nCandFullTops_bdt2->Fill(nBDTTops);
	      break;
	    case 3: 
	      h_nCandFullTops_hot3->Fill(nHOTTops);
	      h_nCandFullTops_bdt3->Fill(nBDTTops);
	      break;
	    case 4: 
	      h_nCandFullTops_hot4->Fill(nHOTTops);
	      h_nCandFullTops_bdt4->Fill(nBDTTops);
	      break;
	    }
            printf("\n");
	}

    //debug counter:
    std::cout << "\n\n\n\n=====================================\nnERROR: " << nERROR << "\n=====================================" << std::endl;

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
    //of = new TFile("results.root", "RECREATE");
    of = new TFile(outFile.c_str(), "RECREATE");

    //Write nJet Histo
    h_nJet->Write();
    h_nJet_N1_hot->Write();
    h_nJet_N2_hot->Write();
    h_nJet_N3_hot->Write();
    h_nJet_N1_bdt->Write();
    h_nJet_N2_bdt->Write();
    h_nJet_N3_bdt->Write();
    
    //Write Event Histos
    h_nTrueRecoTops_full->Write();
    h_nTrueRecoTops->Write();
    h_nCandTops_hot->Write();
    h_nCandTops_bdt->Write();
    h_nCandTops_hot1->Write();
    h_nCandTops_bdt1->Write();
    h_nCandTops_hot2->Write();
    h_nCandTops_bdt2->Write();
    h_nCandTops_hot3->Write();
    h_nCandTops_bdt3->Write();
    h_nCandTops_hot4->Write();
    h_nCandTops_bdt4->Write();
    h_nCandTops_hot_oth->Write();
    h_nCandTops_bdt_oth->Write();

    h_nCandFullTops_hot1->Write();
    h_nCandFullTops_bdt1->Write();
    h_nCandFullTops_hot2->Write();
    h_nCandFullTops_bdt2->Write();
    h_nCandFullTops_hot3->Write();
    h_nCandFullTops_bdt3->Write();
    h_nCandFullTops_hot4->Write();
    h_nCandFullTops_bdt4->Write();
    h_nCandFullTops_hot_0->Write();
    h_nCandFullTops_bdt_0->Write();

    //Write HOT histos
    h_typeIII_hot->Write();
    h_typeIIb_hot->Write();
    h_typeIIw_hot->Write();
    h_typeIImib_hot->Write();
    h_typeIImiq_hot->Write();
    h_typeIt_hot->Write();
    h_typeIp_hot->Write();
    h_type0xI_hot->Write();

    h_eventN1_sum_hot->Write();
    h_eventN1_III_hot->Write();
    h_eventN1_II_hot->Write();
    h_eventN1_IIb_hot->Write();
    h_eventN1_IImib_hot->Write();
    h_eventN1_IIw_hot->Write();
    h_eventN1_IImiq_hot->Write();
    h_eventN1_Itp_hot->Write();
    h_eventN1_It_hot->Write();
    h_eventN1_Ip_hot->Write();
    h_eventN1_0xI_hot->Write();
    h_eventN1_I_hot->Write();
    h_eventN1_0_hot->Write();

    h_eventN2_sum_hot->Write();
    h_eventN2_III_hot->Write();
    h_eventN2_II_hot->Write();
    h_eventN2_IIb_hot->Write();
    h_eventN2_IImib_hot->Write();
    h_eventN2_IIw_hot->Write();
    h_eventN2_IImiq_hot->Write();
    h_eventN2_Itp_hot->Write();
    h_eventN2_It_hot->Write();
    h_eventN2_Ip_hot->Write();
    h_eventN2_0xI_hot->Write();
    h_eventN2_I_hot->Write();
    h_eventN2_0_hot->Write();

    h_eventN3_sum_hot->Write();
    h_eventN3_III_hot->Write();
    h_eventN3_II_hot->Write();
    h_eventN3_IIb_hot->Write();
    h_eventN3_IImib_hot->Write();
    h_eventN3_IIw_hot->Write();
    h_eventN3_IImiq_hot->Write();
    h_eventN3_Itp_hot->Write();
    h_eventN3_It_hot->Write();
    h_eventN3_Ip_hot->Write();
    h_eventN3_0xI_hot->Write();
    h_eventN3_I_hot->Write();
    h_eventN3_0_hot->Write();

    //Write BDT histos
    h_typeIII_bdt->Write();
    h_typeIIb_bdt->Write();
    h_typeIIw_bdt->Write();
    h_typeIImib_bdt->Write();
    h_typeIImiq_bdt->Write();
    h_typeIt_bdt->Write();
    h_typeIp_bdt->Write();
    h_type0xI_bdt->Write();

    h_eventN1_sum_bdt->Write();
    h_eventN1_III_bdt->Write();
    h_eventN1_II_bdt->Write();
    h_eventN1_IIb_bdt->Write();
    h_eventN1_IImib_bdt->Write();
    h_eventN1_IIw_bdt->Write();
    h_eventN1_IImiq_bdt->Write();
    h_eventN1_Itp_bdt->Write();
    h_eventN1_It_bdt->Write();
    h_eventN1_Ip_bdt->Write();
    h_eventN1_0xI_bdt->Write();
    h_eventN1_I_bdt->Write();
    h_eventN1_0_bdt->Write();

    h_eventN2_sum_bdt->Write();
    h_eventN2_III_bdt->Write();
    h_eventN2_II_bdt->Write();
    h_eventN2_IIb_bdt->Write();
    h_eventN2_IImib_bdt->Write();
    h_eventN2_IIw_bdt->Write();
    h_eventN2_IImiq_bdt->Write();
    h_eventN2_Itp_bdt->Write();
    h_eventN2_It_bdt->Write();
    h_eventN2_Ip_bdt->Write();
    h_eventN2_0xI_bdt->Write();
    h_eventN2_I_bdt->Write();
    h_eventN2_0_bdt->Write();

    h_eventN3_sum_bdt->Write();
    h_eventN3_III_bdt->Write();
    h_eventN3_II_bdt->Write();
    h_eventN3_IIb_bdt->Write();
    h_eventN3_IImib_bdt->Write();
    h_eventN3_IIw_bdt->Write();
    h_eventN3_IImiq_bdt->Write();
    h_eventN3_Itp_bdt->Write();
    h_eventN3_It_bdt->Write();
    h_eventN3_Ip_bdt->Write();
    h_eventN3_0xI_bdt->Write();
    h_eventN3_I_bdt->Write();
    h_eventN3_0_bdt->Write();
    
    //Special top-jet counting histograms

    h_eventN1_CntJet_hot->Write();  
    h_eventN2_CntJet_hot->Write();  
    h_eventN3_CntJet_hot->Write();  
    h_eventN12_CntJet_hot->Write(); 
    h_eventN123_CntJet_hot->Write();
                           
    h_eventN1_CntJet_bdt->Write();  
    h_eventN2_CntJet_bdt->Write();  
    h_eventN3_CntJet_bdt->Write();  
    h_eventN12_CntJet_bdt->Write(); 
    h_eventN123_CntJet_bdt->Write();

    //Special BDT cross-check histograms
    h_noscaling_bdt->Write();           
    h_noscaling2_bdt->Write();
    h_scaling_bdt->Write();
    h_scaling2_bdt->Write();

    h_nERROR->Write();

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

    // Event histos
    delete h_nJet;
    delete h_nJet_N1_hot;
    delete h_nJet_N2_hot;
    delete h_nJet_N3_hot;
    delete h_nJet_N1_bdt;
    delete h_nJet_N2_bdt;
    delete h_nJet_N3_bdt;

    delete h_nTrueRecoTops_full;
    delete h_nTrueRecoTops;
    delete h_nCandTops_hot;	       
    delete h_nCandTops_bdt;
    delete h_nCandTops_hot1;
    delete h_nCandTops_bdt1;
    delete h_nCandTops_hot2;
    delete h_nCandTops_bdt2;
    delete h_nCandTops_hot3;
    delete h_nCandTops_bdt3;
    delete h_nCandTops_hot4;
    delete h_nCandTops_bdt4;
    delete h_nCandTops_hot_oth;
    delete h_nCandTops_bdt_oth;

    delete h_nCandFullTops_hot1;
    delete h_nCandFullTops_bdt1;
    delete h_nCandFullTops_hot2;
    delete h_nCandFullTops_bdt2;
    delete h_nCandFullTops_hot3;
    delete h_nCandFullTops_bdt3;
    delete h_nCandFullTops_hot4;
    delete h_nCandFullTops_bdt4;
    delete h_nCandFullTops_hot_0;
    delete h_nCandFullTops_bdt_0;

    // HOT histos
    delete h_typeIII_hot;
    delete h_typeIIb_hot;
    delete h_typeIIw_hot;
    delete h_typeIImib_hot;
    delete h_typeIImiq_hot;
    delete h_typeIt_hot;
    delete h_typeIp_hot;
    delete h_type0xI_hot;

    delete h_eventN1_sum_hot;
    delete h_eventN1_III_hot;
    delete h_eventN1_II_hot;
    delete h_eventN1_Itp_hot;
    delete h_eventN1_0xI_hot;
    delete h_eventN1_IIb_hot;
    delete h_eventN1_IImib_hot;
    delete h_eventN1_IIw_hot;
    delete h_eventN1_IImiq_hot;
    delete h_eventN1_It_hot;
    delete h_eventN1_Ip_hot;
    delete h_eventN1_I_hot;
    delete h_eventN1_0_hot;

    delete h_eventN2_sum_hot;
    delete h_eventN2_III_hot;
    delete h_eventN2_II_hot;
    delete h_eventN2_Itp_hot;
    delete h_eventN2_0xI_hot;
    delete h_eventN2_IIb_hot;
    delete h_eventN2_IImib_hot;
    delete h_eventN2_IIw_hot;
    delete h_eventN2_IImiq_hot;
    delete h_eventN2_It_hot;
    delete h_eventN2_Ip_hot;
    delete h_eventN2_I_hot;
    delete h_eventN2_0_hot;

    delete h_eventN3_sum_hot;
    delete h_eventN3_III_hot;
    delete h_eventN3_II_hot;
    delete h_eventN3_Itp_hot;
    delete h_eventN3_0xI_hot;
    delete h_eventN3_IIb_hot;
    delete h_eventN3_IImib_hot;
    delete h_eventN3_IIw_hot;
    delete h_eventN3_IImiq_hot;
    delete h_eventN3_It_hot;
    delete h_eventN3_Ip_hot;
    delete h_eventN3_I_hot;
    delete h_eventN3_0_hot;


    // BDT histos
    delete h_typeIII_bdt;
    delete h_typeIIb_bdt;
    delete h_typeIIw_bdt;
    delete h_typeIImib_bdt;
    delete h_typeIImiq_bdt;
    delete h_typeIt_bdt;
    delete h_typeIp_bdt;
    delete h_type0xI_bdt;

    delete h_eventN1_sum_bdt;
    delete h_eventN1_III_bdt;
    delete h_eventN1_II_bdt;
    delete h_eventN1_Itp_bdt;
    delete h_eventN1_0xI_bdt;
    delete h_eventN1_IIb_bdt;
    delete h_eventN1_IImib_bdt;
    delete h_eventN1_IIw_bdt;
    delete h_eventN1_IImiq_bdt;
    delete h_eventN1_It_bdt;
    delete h_eventN1_Ip_bdt;
    delete h_eventN1_I_bdt;
    delete h_eventN1_0_bdt;

    delete h_eventN2_sum_bdt;
    delete h_eventN2_III_bdt;
    delete h_eventN2_II_bdt;
    delete h_eventN2_Itp_bdt;
    delete h_eventN2_0xI_bdt;
    delete h_eventN2_IIb_bdt;
    delete h_eventN2_IImib_bdt;
    delete h_eventN2_IIw_bdt;
    delete h_eventN2_IImiq_bdt;
    delete h_eventN2_It_bdt;
    delete h_eventN2_Ip_bdt;
    delete h_eventN2_I_bdt;
    delete h_eventN2_0_bdt;

    delete h_eventN3_sum_bdt;
    delete h_eventN3_III_bdt;
    delete h_eventN3_II_bdt;
    delete h_eventN3_Itp_bdt;
    delete h_eventN3_0xI_bdt;
    delete h_eventN3_IIb_bdt;
    delete h_eventN3_IImib_bdt;
    delete h_eventN3_IIw_bdt;
    delete h_eventN3_IImiq_bdt;
    delete h_eventN3_It_bdt;
    delete h_eventN3_Ip_bdt;
    delete h_eventN3_I_bdt;
    delete h_eventN3_0_bdt;

    delete h_noscaling_bdt;
    delete h_noscaling2_bdt;
    delete h_scaling_bdt;
    delete h_scaling2_bdt;

    delete h_eventN1_CntJet_hot;
    delete h_eventN2_CntJet_hot;
    delete h_eventN3_CntJet_hot;
    delete h_eventN12_CntJet_hot;
    delete h_eventN123_CntJet_hot;
                       
    delete h_eventN1_CntJet_bdt;
    delete h_eventN2_CntJet_bdt;
    delete h_eventN3_CntJet_bdt;
    delete h_eventN12_CntJet_bdt;
    delete h_eventN123_CntJet_bdt;

    delete h_nERROR;
    delete tf, tf2, td, tree, of;
    //delete h_typeIII, h_typeIIb, h_typeIIw, h_typeIIo, h_typeI, RecoTypes;
    exit(0);
}
