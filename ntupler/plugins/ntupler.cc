// -*- C++ -*-
//
// Package:    Demo/ntupler
// Class:      ntupler
// 
/**\class ntupler ntupler.cc Demo/ntupler/plugins/ntupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nick Manganelli
//         Created:  Tue, 28 Aug 2018 15:30:59 GMT
//
//

//Generate dictionary for vector<TLorentsVector> so that ROOT will be happy making/streaming such a branch
//#ifdef __CINT__
//#pragma link C++ class std::vector<TLorentzVector>+;
//#endif
//Didn't work...


// system include files
#include <memory> //e.g. for unique_ptr made with make_unique
#include <cmath>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
//#include <Rtypes.h>
//#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h" //threadsafe option
//#include "FWCore/Framework/interface/EDAnalyzer.h" //thread UNSAFE option (TFileService?)

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class ntupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
  //explicit 
      explicit ntupler(const edm::ParameterSet&);
      ~ntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  bool isData, isMC;
  bool is2016, is2017, is2018;

  //Tokens
  edm::EDGetTokenT<std::vector<pat::Jet>> JetToken;
  edm::EDGetTokenT<std::vector<pat::Muon>> MuonToken;
  edm::EDGetTokenT<std::vector<pat::Electron>> ElectronToken;
  edm::EDGetTokenT<std::vector<pat::MET>> METToken;
  edm::EDGetTokenT<std::vector<reco::Vertex>> VtxToken;
  edm::EDGetTokenT<edm::TriggerResults> HLTToken;
  edm::EDGetTokenT<edm::TriggerResults> FltToken;

  //HLT Triggers
  std::vector<std::string> HLT_MuMu_S;
  std::vector<std::string> HLT_ElMu_S;
  std::vector<std::string> HLT_ElEl_S;
  std::vector<std::string> HLT_Mu_S;
  std::vector<std::string> HLT_El_S;

  //Use deque instead of vector, which for type bool has a specialized format (8 bools per byte) and whose elements don't behave as C++ bools (one BYTE each)
  boost::dynamic_bitset<> HLT_MuMu_B;
  boost::dynamic_bitset<> HLT_ElMu_B;
  boost::dynamic_bitset<> HLT_ElEl_B;
  boost::dynamic_bitset<> HLT_Mu_B;
  boost::dynamic_bitset<> HLT_El_B;

  //TTree
  TTree *tree;

  // TBranch *br_my_var;
  // TBranch *br_nEvts;
  // TBranch *br_my_other_var;


  uint nEvts, run, lumiBlock, nEvent;
  //bool HLT
  bool MuMu, ElMu, ElEl, El, Mu, SL, DL;
  std::vector<TLorentzVector> *JetVec, *MuonVec, *ElectronVec;
  std::vector<double> *qgPtDVec, *qgAxis1Vec, *qgAxis2Vec, *qgMultVec;
  std::vector<double> *deepCSVbVec, *deepCSVcVec, *deepCSVlVec, *deepCSVbbVec, *deepCSVccVec, *btagVec;
  std::vector<double> *chargedHadronEnergyFractionVec, *neutralHadronEnergyFractionVec, *chargedEmEnergyFractionVec;
  std::vector<double> *neutralEmEnergyFractionVec, *muonEnergyFractionVec, *photonEnergyFractionVec, *electronEnergyFractionVec;
  std::vector<double> *recoJetsHFHadronEnergyFractionVec, *recoJetsHFEMEnergyFractionVec;
  std::vector<double> *chargedHadronMultiplicityVec, *neutralHadronMultiplicityVec, *photonMultiplicityVec, *electronMultiplicityVec, *muonMultiplicityVec;

  TBranch *br_nEvts2, *br_run, *br_lumiBlock, *br_nEvent;
  TBranch *br_JetVec, *br_MuonVec, *br_ElectronVec;

  //Jet sub-branches for HOT tagger
  TBranch *br_qgPtD, *br_qgAxis1, *br_qgAxis2, *br_qgMult;
  TBranch *br_deepCSVb, *br_deepCSVc, *br_deepCSVl, *br_deepCSVbb, *br_deepCSVcc, *br_btag;
  TBranch *br_chargedHadronEnergyFraction, *br_neutralHadronEnergyFraction, *br_chargedEmEnergyFraction;
  TBranch *br_neutralEmEnergyFraction, *br_muonEnergyFraction, *br_photonEnergyFraction, *br_electronEnergyFraction;
  TBranch *br_recoJetsHFHadronEnergyFraction, *br_recoJetsHFEMEnergyFraction;
  TBranch *br_chargedHadronMultiplicity, *br_neutralHadronMultiplicity, *br_photonMultiplicity, *br_electronMultiplicity, *br_muonMultiplicity;
};

//
// constants, enums and typedefs
//
  // MET Filters
  //std::string HBHENoiseFilter_Selector_;
  //std::string EEBadScNoiseFilter_Selector_;
//
// static data member definitions
//

//
// constructors and destructor
//
ntupler::ntupler(const edm::ParameterSet& iConfig)//:nEvts(0)//, my_var(0)
{
   //Explicitly declare shared resource TFileService to make it threadsafe
   usesResource("TFileService");

   //FIXME : HardCoding for bools
   isData = false;
   isMC = true;
   is2016 = true;
   is2017 = false;
   is2018 = false;


   JetToken = consumes<std::vector<pat::Jet> >(edm::InputTag("selectedUpdatedPatJetsDeepCSV"));
   MuonToken = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
   ElectronToken = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons"));
   METToken = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));
   VtxToken = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices")); //FIXME
   HLTToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
   FltToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO"));

   //MET Filter Settings
   //HBHENoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("HBHENoiseFilter_Selector_");
   //EEBadScNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("EEBadScNoiseFilter_Selector_");
   //HBHENoiseFilter_Selector_ = "HBHENoiseFilter_Selector_";
   //EEBadScNoiseFilter_Selector_ = "EEBadScNoiseFilter_Selector_";

   if(is2016){
     std::cout << "Defining triggers for 2016" << std::endl;
     //MuMu Triggers
     HLT_MuMu_S.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6");
     HLT_MuMu_S.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7");
     HLT_MuMu_S.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5");
     HLT_MuMu_S.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6");

     //ElMu Triggers
     HLT_ElMu_S.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v7");
     HLT_ElMu_S.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v3");
     HLT_ElMu_S.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9");
     HLT_ElMu_S.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4");

     //ElEl Triggers
     HLT_ElEl_S.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9");

     //Mu Triggers
     HLT_Mu_S.push_back("HLT_IsoTkMu24_v4");
     HLT_Mu_S.push_back("HLT_IsoMu24_v4");

     //El Triggers
     HLT_El_S.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
   }
   else if(is2017){
     std::cout << "FIX ME!" << std::endl;
   }
   else if(is2018){
     std::cout << "FIX ME!" << std::endl;
   }
   else
     std::cout << "Error: Data is not from 2016, 2017, or 2018. Is it ReReco?" << std::endl;

   //Set HLT bits to zero via initialization, using the boost::dynamic_bitset, with the proper number of bits
   HLT_MuMu_B = boost::dynamic_bitset<>(HLT_MuMu_S.size(), 0ul);
   HLT_ElMu_B = boost::dynamic_bitset<>(HLT_ElMu_S.size(), 0ul);
   HLT_ElEl_B = boost::dynamic_bitset<>(HLT_ElEl_S.size(), 0ul);
   HLT_Mu_B = boost::dynamic_bitset<>(HLT_Mu_S.size(), 0ul);
   HLT_El_B = boost::dynamic_bitset<>(HLT_El_S.size(), 0ul);

}


ntupler::~ntupler()
{
   usesResource("TFileService");
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<std::vector<pat::Jet> > jets;
   iEvent.getByToken(JetToken, jets);
   edm::Handle<std::vector<pat::Muon> > muons;
   iEvent.getByToken(MuonToken, muons);
   edm::Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByToken(ElectronToken, electrons);
   edm::Handle<std::vector<pat::MET> > mets;
   iEvent.getByToken(METToken, mets);
   edm::Handle<edm::TriggerResults> METFlt;
   iEvent.getByToken(FltToken, METFlt);
   edm::Handle<edm::TriggerResults> HLTTrg;
   iEvent.getByToken(HLTToken, HLTTrg);

   //See /afs/cern.ch/user/n/nmangane/DAS/EGammaExercises/MuonExercise3/plugins/MuonExercise3.cc for histogram array creation, tight (vertex) selection
   //edm::Handle <pat::PackedGenParticleCollection> genColl;
   //iEvent.getByToken(genCollToken, genColl);
   //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ CP code
  // // Retrieve the GenParticle collection from the event 
  // edm::Handle<reco::GenParticleCollection> genParticles;
  // iEvent.getByToken(genPartToken, genParticles);
  // std::vector<reco::GenParticle> genColl = (*genParticles);

  // // Just for peace of mind, let's check that the collection is valid... 
  // if(!genParticles.isValid()) {
  //   throw cms::Exception("GenParticle collection not valid!"); 
  // } 

  // // Same with pat::Muons 
  // edm::Handle<edm::View<pat::Muon> > muons;
  // iEvent.getByToken(muonToken, muons);
  // if(!muons.isValid()) {
  //   throw cms::Exception("Muon collection not valid!"); 
  // }

  // // Same with vertices 
  // edm::Handle<std::vector<reco::Vertex>> vertices;
  // iEvent.getByToken(vertexToken, vertices);
  // if(!vertices.isValid()) {
  //   throw cms::Exception("Vertex collection not valid!"); 
  // } 

  // // Let's check that we have at least one good vertex! 
  // std::vector<reco::Vertex>::const_iterator firstGoodVertex = vertices->end();

  // for (std::vector<reco::Vertex>::const_iterator it=vertices->begin(); it!=firstGoodVertex; ++it) {
  //   if (!it->isFake() && it->ndof()>4 && it->position().Rho()<2. && std::abs(it->position().Z())<24.) {
  //     if(firstGoodVertex == vertices->end()) firstGoodVertex = it;
  //     break;
  //   }
  // }

  // // Require a good vertex
  // if(firstGoodVertex == vertices->end()) return;
  // // int nvtx = vertices->size(); 

  // //std::cout << deb++ << std::endl;
  // edm::View<pat::Muon>::const_iterator muend = muons->end();
  // for (edm::View<pat::Muon>::const_iterator it=muons->begin(); it!=muend; ++it) {
  //   // Require muon to have a silicon track (i.e. global- or tracker-muon) 
  //   if (!it->innerTrack().isNonnull()) continue;

  //   // Require that muon be within the tracker volume and have pt > 20 GeV
  //   if (it->pt() < 20 || it->eta() > 2.4) continue;
    
  //   // Let's check the origin of the muon: prompt, HF, LF, other?
  //   MuonParentage parentage = MuonParentage::OTHER;

  //   // For this, we need to find the gen particle associated with our reconstructed muon
  //   // Since we are using pat::Muons, the gen-matching is already done for us!
  //   // No need to loop over the GenParticle collection and perform a geometrical matching 
  //   const reco::GenParticle* gp = it->genParticle();

  //   // Check if the pat::Muon has a GenParticle associated
  //   // In what cases there is no gen-matching? 
  //   if(gp!=0) {
  //     // This function determines the muon origin for you
  //     // Take some time to understand what it does 
  //     parentage = getParentType(*gp); 
  //   }
  //   else {
  //     // If there is no genParticle() and it's a Drell-Yan or top-top sample, stop here!
  //     // Proceed with the classification only when running on the QCD sample! 
  //     // In all the other samples, muons from light flavor decays are NOT saved
  //     // in the GenParticle collection. Therefore, light-flavor decays would be
  //     // classified as "other", which is not correct.
  //     // QCD samples, on the other hand, have gen-particles also for light-flavor decays 
  //     if(!isQCD) continue; 
  //   }
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ end CP code

   // HLT Triggers
   // WARNING: version will differ in different productions, must be parsed in the future to make this more automated.
   //Store the string name of triggers in arrays with postfix _S, and corresponding booleans will be in arrays with postfix _B


   const edm::TriggerNames &HLTnames = iEvent.triggerNames(*HLTTrg);
   for (unsigned int i = 0, n = HLTTrg->size(); i < n; ++i) {
     //std::cout << HLTnames.triggerName(i) << std::endl;
     //FIXME Option: Create a mapping on a per-year basis to try matching based on position, and if that fails, then loop through instead
     //For efficicency: Do one search through names to find locations, then store them in global variables to be used in all following events

     for(uint j = 0; j < HLT_MuMu_S.size(); j++)
       if (HLTnames.triggerName(i) == HLT_MuMu_S[j]) HLT_MuMu_B[j] = HLTTrg->accept(i);
     for(uint jj = 0; jj < HLT_ElMu_S.size(); jj++)
       if (HLTnames.triggerName(i) == HLT_ElMu_S[jj]) HLT_ElMu_B[jj] = HLTTrg->accept(i);
     for(uint jjj = 0; jjj < HLT_ElEl_S.size(); jjj++)
       if (HLTnames.triggerName(i) == HLT_ElEl_S[jjj]) HLT_ElEl_B[jjj] = HLTTrg->accept(i);
     for(uint k = 0; k < HLT_Mu_S.size(); k++)
       if (HLTnames.triggerName(i) == HLT_Mu_S[k]) HLT_Mu_B[k] = HLTTrg->accept(i);
     for(uint kk = 0; kk < HLT_ElEl_S.size(); kk++)
       if (HLTnames.triggerName(i) == HLT_El_S[kk]) HLT_El_B[kk] = HLTTrg->accept(i);
     std::cout << HLTnames.triggerName(i) << " " << HLT_MuMu_B << " " << HLTTrg->accept(i) << std::endl;

     // //testing
     // std::cout << HLT_MuMu_B[0] << " " << HLT_MuMu_B[1] << " " << HLT_MuMu_B[2] << " " << HLT_MuMu_B[3] << " " << HLT_ElMu_B[0] << " " << HLT_ElMu_B[1]
     // 	       << " " << HLT_ElMu_B[2] << " " << HLT_ElMu_B[3] << " " << HLT_ElEl_B[0] << " " << HLT_Mu_B[0] << " " << HLT_Mu_B[1] << " " << HLT_El_B[0] << std::endl;
     
   }

   // MET Filters
   const pat::MET &met = mets->front();
   //int passFilterHBHE=0;
   //int passFilterEEBadSC=0;
   const edm::TriggerNames &names = iEvent.triggerNames(*METFlt);
   for (unsigned int i = 0, n = METFlt->size(); i < n; ++i) {
     std::cout << names.triggerName(i) << std::endl;
     //if (names.triggerName(i) == HBHENoiseFilter_Selector_)
     //  passFilterHBHE=metflt->accept(i);
     //if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
     //  passFilterEEBadSC=metflt->accept(i);
   }
   // End filters stuff


   ////////////////////
   //// Event info ////
   ////////////////////
   run = iEvent.id().run();
   lumiBlock = iEvent.id().luminosityBlock();
   nEvent = iEvent.id().event();

   ////////////////////////
   //// Selected Muons ////
   ////////////////////////
   for(const pat::Muon& muon : *muons){
     TLorentzVector perMuonLVec;
     perMuonLVec.SetPtEtaPhiE( muon.pt(), muon.eta(), muon.phi(), muon.energy() );
     MuonVec->push_back(perMuonLVec);
   }

   ///////////////////////
   //// Selected Jets ////
   ///////////////////////
   for(const pat::Jet&jet : *jets){
     TLorentzVector perJetLVec;
     perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
     JetVec->push_back(perJetLVec);

     double qgPtD = jet.userFloat("QGTagger:ptD");
     double qgAxis1 = jet.userFloat("QGTagger:axis1");
     double qgAxis2 = jet.userFloat("QGTagger:axis2");
     double qgMult = static_cast<double>(jet.userInt("QGTagger:mult"));
     double deepCSVb = jet.bDiscriminator("pfDeepCSVJetTags:probb");
     double deepCSVc = jet.bDiscriminator("pfDeepCSVJetTags:probc");
     double deepCSVl = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
     double deepCSVbb = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
     double deepCSVcc = jet.bDiscriminator("pfDeepCSVJetTags:probcc");
     double btag = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
     double chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
     double neutralHadronEnergyFraction = jet.neutralHadronEnergyFraction();
     double chargedEmEnergyFraction = jet.chargedEmEnergyFraction();
     double neutralEmEnergyFraction = jet.neutralEmEnergyFraction();
     double muonEnergyFraction = jet.muonEnergyFraction();
     double photonEnergyFraction = jet.photonEnergyFraction();
     double electronEnergyFraction = jet.electronEnergyFraction();
     double recoJetsHFHadronEnergyFraction = jet.HFHadronEnergyFraction();
     double recoJetsHFEMEnergyFraction = jet.HFEMEnergyFraction();
     double chargedHadronMultiplicity = jet.chargedHadronMultiplicity();
     double neutralHadronMultiplicity = jet.neutralHadronMultiplicity();
     double photonMultiplicity = jet.photonMultiplicity();
     double electronMultiplicity = jet.electronMultiplicity();
     double muonMultiplicity = jet.muonMultiplicity();
    
     qgPtDVec->push_back(qgPtD);
     qgAxis1Vec->push_back(qgAxis1);
     qgAxis2Vec->push_back(qgAxis2);
     qgMultVec->push_back(qgMult);
     deepCSVbVec->push_back(deepCSVb);
     deepCSVcVec->push_back(deepCSVc);
     deepCSVlVec->push_back(deepCSVl);
     deepCSVbbVec->push_back(deepCSVbb);
     deepCSVccVec->push_back(deepCSVcc);
     btagVec->push_back(btag);
     chargedHadronEnergyFractionVec->push_back(chargedHadronEnergyFraction);
     neutralHadronEnergyFractionVec->push_back(neutralHadronEnergyFraction);
     chargedEmEnergyFractionVec->push_back(chargedEmEnergyFraction);
     neutralEmEnergyFractionVec->push_back(neutralEmEnergyFraction);
     muonEnergyFractionVec->push_back(muonEnergyFraction);
     photonEnergyFractionVec->push_back(photonEnergyFraction);
     electronEnergyFractionVec->push_back(electronEnergyFraction);
     recoJetsHFHadronEnergyFractionVec->push_back(recoJetsHFHadronEnergyFraction);
     recoJetsHFEMEnergyFractionVec->push_back(recoJetsHFEMEnergyFraction);
     chargedHadronMultiplicityVec->push_back(chargedHadronMultiplicity);
     neutralHadronMultiplicityVec->push_back(neutralHadronMultiplicity);
     photonMultiplicityVec->push_back(photonMultiplicity);
     electronMultiplicityVec->push_back(electronMultiplicity);
     muonMultiplicityVec->push_back(muonMultiplicity);
 
     std::cout << " Pt: " << jet.pt() << " btag: " << btag << std::endl;
   }

   nEvts++;
   tree->Fill();

   //Clear pointers
   JetVec->clear();
   MuonVec->clear();
   ElectronVec->clear();
   qgPtDVec->clear();
   qgAxis1Vec->clear();
   qgAxis2Vec->clear();
   qgMultVec->clear();
   deepCSVbVec->clear();
   deepCSVcVec->clear();
   deepCSVlVec->clear();
   deepCSVbbVec->clear();
   deepCSVccVec->clear();
   btagVec->clear();
   chargedHadronEnergyFractionVec->clear();
   neutralHadronEnergyFractionVec->clear();
   chargedEmEnergyFractionVec->clear();
   neutralEmEnergyFractionVec->clear();
   muonEnergyFractionVec->clear();
   photonEnergyFractionVec->clear();
   electronEnergyFractionVec->clear();
   recoJetsHFHadronEnergyFractionVec->clear();
   recoJetsHFEMEnergyFractionVec->clear();
   chargedHadronMultiplicityVec->clear();
   neutralHadronMultiplicityVec->clear();
   photonMultiplicityVec->clear();
   electronMultiplicityVec->clear();
   muonMultiplicityVec->clear();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ntupler::beginJob()
{   
   nEvts = 0;
   run = -1;
   lumiBlock = -1; 
   nEvent = -1;
   // HLT__= -1;
   // HLT__= -1;
   // HLT__= -1;
   // HLT__= -1;
   // HLT__= -1;
   //FIXME
   JetVec = new std::vector<TLorentzVector>;
   MuonVec = new std::vector<TLorentzVector>;
   ElectronVec = new std::vector<TLorentzVector>;
   qgPtDVec = new std::vector<double>;
   qgAxis1Vec = new std::vector<double>;
   qgAxis2Vec = new std::vector<double>;
   qgMultVec = new std::vector<double>;
   deepCSVbVec = new std::vector<double>;
   deepCSVcVec = new std::vector<double>;
   deepCSVlVec = new std::vector<double>;
   deepCSVbbVec = new std::vector<double>;
   deepCSVccVec = new std::vector<double>;
   btagVec = new std::vector<double>;
   chargedHadronEnergyFractionVec = new std::vector<double>;
   neutralHadronEnergyFractionVec = new std::vector<double>;
   chargedEmEnergyFractionVec = new std::vector<double>;
   neutralEmEnergyFractionVec = new std::vector<double>;
   muonEnergyFractionVec = new std::vector<double>;
   photonEnergyFractionVec = new std::vector<double>;
   electronEnergyFractionVec = new std::vector<double>;
   recoJetsHFHadronEnergyFractionVec = new std::vector<double>;
   recoJetsHFEMEnergyFractionVec = new std::vector<double>;
   chargedHadronMultiplicityVec = new std::vector<double>;
   neutralHadronMultiplicityVec = new std::vector<double>;
   photonMultiplicityVec = new std::vector<double>;
   electronMultiplicityVec = new std::vector<double>;
   muonMultiplicityVec = new std::vector<double>;

   usesResource("TFileService");
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("nTuple"," Event nTuple");
   br_run = tree->Branch("run", &run);
   br_lumiBlock = tree->Branch("lumiBlock", &lumiBlock);
   br_nEvent = tree->Branch("event", &nEvent);
   br_JetVec = tree->Branch("JetVec", "vector<TLorentzVector>", &JetVec, 32000,-1);
   br_MuonVec = tree->Branch("MuonVec", "vector<TLorentzVector>", &MuonVec, 32000,-1);
   br_ElectronVec = tree->Branch("ElectronVec", "vector<TLorentzVector>", &ElectronVec, 32000,-1);
   tree->Branch("qgPtD", &qgPtDVec);
   tree->Branch("qgAxis1", &qgAxis1Vec);
   tree->Branch("qgAxis2", &qgAxis2Vec);
   tree->Branch("qgMult", &qgMultVec);
   tree->Branch("deepCSVb", &deepCSVbVec); 
   tree->Branch("deepCSVc", &deepCSVcVec); 
   tree->Branch("deepCSVl", &deepCSVlVec);
   tree->Branch("deepCSVbb", &deepCSVbbVec);
   tree->Branch("deepCSVcc", &deepCSVccVec);
   tree->Branch("btag", &btagVec);
   tree->Branch("chargedHadronEnergyFraction", &chargedHadronEnergyFractionVec);
   tree->Branch("neutralHadronEnergyFraction", &neutralHadronEnergyFractionVec);
   tree->Branch("chargedEmEnergyFraction", &chargedEmEnergyFractionVec);
   tree->Branch("neutralEmEnergyFraction", &neutralEmEnergyFractionVec);
   tree->Branch("muonEnergyFraction", &muonEnergyFractionVec);
   tree->Branch("photonEnergyFraction", &photonEnergyFractionVec);
   tree->Branch("electronEnergyFraction", &electronEnergyFractionVec);
   tree->Branch("recoJetsHFHadronEnergyFraction", &recoJetsHFHadronEnergyFractionVec);
   tree->Branch("recoJetsHFEMEnergyFraction", &recoJetsHFEMEnergyFractionVec);
   tree->Branch("chargedHadronMultiplicity", &chargedHadronMultiplicityVec);
   tree->Branch("neutralHadronMultiplicity", &neutralHadronMultiplicityVec);
   tree->Branch("photonMultiplicity", &photonMultiplicityVec);
   tree->Branch("electronMultiplicity", &electronMultiplicityVec);
   tree->Branch("muonMultiplicity", &muonMultiplicityVec);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntupler::endJob() 
{
   usesResource("TFileService");
   //********** WARNING **********//
   //Using tree->Write(***) will invoke Write() TWICE, as the cmsRun framework does this for you. This results in a duplicate tree in the file
   // tree->GetDirectory()->cd();
   // tree->Write("", TObject::kOverwrite);
   delete JetVec;
   delete MuonVec;
   delete ElectronVec;
   delete qgPtDVec;
   delete qgAxis1Vec;
   delete qgAxis2Vec;
   delete qgMultVec;
   delete deepCSVbVec; 
   delete deepCSVcVec; 
   delete deepCSVlVec;
   delete deepCSVbbVec;
   delete deepCSVccVec;
   delete btagVec;
   delete chargedHadronEnergyFractionVec;
   delete neutralHadronEnergyFractionVec;
   delete chargedEmEnergyFractionVec;
   delete neutralEmEnergyFractionVec; 
   delete muonEnergyFractionVec; 
   delete photonEnergyFractionVec; 
   delete electronEnergyFractionVec;
   delete recoJetsHFHadronEnergyFractionVec; 
   delete recoJetsHFEMEnergyFractionVec;
   delete chargedHadronMultiplicityVec; 
   delete neutralHadronMultiplicityVec; 
   delete photonMultiplicityVec; 
   delete electronMultiplicityVec; 
   delete muonMultiplicityVec;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntupler);
