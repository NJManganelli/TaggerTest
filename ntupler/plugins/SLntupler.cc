// -*- C++ -*-
//
// Package:    Demo/SLntupler
// Class:      SLntupler
// 
/**\class SLntupler SLntupler.cc Demo/ntupler/plugins/SLntupler.cc

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
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" //for isolation, not needed
//#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

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

class SLntupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
  //explicit 
      explicit SLntupler(const edm::ParameterSet&);
      ~SLntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  //Internal counter
  uint counter, theProblemEvent;

  //Parameter Set
  bool isData, isMC, deBug, verBose, maskDeepCSV;
  bool is2016, is2017, is2018;
  double HTMin;
  int NjMin;

  //Tokens
  edm::EDGetTokenT<std::vector<pat::Jet> > JetToken;
  edm::EDGetTokenT<std::vector<pat::Muon> > MuonToken;
  edm::EDGetTokenT<std::vector<pat::Electron> > ElectronToken;
  edm::EDGetTokenT<edm::View<reco::GsfElectron> > GsfElectronToken;
  edm::EDGetTokenT<edm::ValueMap<bool> > EleVetoIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool> > EleLooseIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool> > EleMediumIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool> > EleTightIdMapToken;
  edm::EDGetTokenT<std::vector<pat::MET> > METToken;
  edm::EDGetTokenT<std::vector<reco::Vertex> > VtxToken;
  edm::EDGetTokenT<edm::TriggerResults> HLTToken;
  edm::EDGetTokenT<edm::TriggerResults> FltToken;
  edm::EDGetTokenT<std::vector<reco::Conversion> > ConversionsToken;
  edm::EDGetTokenT<reco::BeamSpot> BeamSpotToken;
  edm::EDGetTokenT<double> RhoToken;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> > GenToken;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > PrunedGenToken;
  edm::EDGetTokenT<std::vector<reco::GenJet> > GenJetToken;
  edm::EDGetTokenT<std::vector<reco::GenJet> > GenJetAK8Token;

  //HLT Triggers
  std::vector<std::string> HLT_MuMu_S;
  std::vector<std::string> HLT_ElMu_S;
  std::vector<std::string> HLT_ElEl_S;
  std::vector<std::string> HLT_Mu_S;
  std::vector<std::string> HLT_El_S;

  //Use dynamic bitset instead of vector<bool>, which has a specialized format (8 bools per byte - 1 bit each) and whose elements don't behave as C++ bools (one BYTE each)
  boost::dynamic_bitset<> HLT_MuMu_B;
  boost::dynamic_bitset<> HLT_ElMu_B;
  boost::dynamic_bitset<> HLT_ElEl_B;
  boost::dynamic_bitset<> HLT_Mu_B;
  boost::dynamic_bitset<> HLT_El_B;

  //For storing the bits in ROOT, using explicit casting bitset -> unsigned long, then implicit casting ulong -> uint
  uint HLT_MuMu_Bits, HLT_ElMu_Bits, HLT_ElEl_Bits, HLT_Mu_Bits, HLT_El_Bits;

  //MET Filters
  std::vector<std::string> MET_Flt_S;

  //dynamic bitset for MET Filters
  boost::dynamic_bitset<> MET_Flt_B;

  //for storing filter bits in ROOT
  uint MET_Flt_Bits;

  //TTree
  TTree *tree;

  //Data and collections
  uint nEvts, nRun, nLumiBlock, nEvent;
  bool MuMu, ElMu, ElEl, El, Mu, SL, DL;   //bool HLT
  bool selectedLepIsMu, vetoLep1IsMu, vetoLep2IsMu; //FIXME add these to tree, etc...
  double HT, HTX, HT2M; //FIXME Calculate and add these...
  std::vector<TLorentzVector> *JetLVec, *selectedLepLVec, *vetoLepLVec, *VertexVec, *METLVec; //Change to pointer vector of pointers for sorting efficiency! FIXME!
  std::vector<double> *qgPtDVec, *qgAxis1Vec, *qgAxis2Vec, *qgMultVec;
  std::vector<double> *deepCSVbVec, *deepCSVcVec, *deepCSVlVec, *deepCSVbbVec, *deepCSVccVec, *btagVec;
  std::vector<double> *chargedHadronEnergyFractionVec, *neutralHadronEnergyFractionVec, *chargedEmEnergyFractionVec;
  std::vector<double> *neutralEmEnergyFractionVec, *muonEnergyFractionVec, *photonEnergyFractionVec, *electronEnergyFractionVec;
  std::vector<double> *recoJetsHFHadronEnergyFractionVec, *recoJetsHFEMEnergyFractionVec;
  std::vector<double> *chargedHadronMultiplicityVec, *neutralHadronMultiplicityVec, *photonMultiplicityVec, *electronMultiplicityVec, *muonMultiplicityVec;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SLntupler::SLntupler(const edm::ParameterSet& iConfig): 
  theProblemEvent(iConfig.getParameter<uint>("theProblemEvent")),
  isData(iConfig.getParameter<bool>("isData")), isMC(iConfig.getParameter<bool>("isMC")), 
  deBug(iConfig.getParameter<bool>("deBug")), verBose(iConfig.getParameter<bool>("verBose")),
  maskDeepCSV(iConfig.getParameter<bool>("maskDeepCSV")), is2016(iConfig.getParameter<bool>("is2016")), 
  is2017(iConfig.getParameter<bool>("is2017")), is2018(iConfig.getParameter<bool>("is2018")), 
  HTMin(iConfig.getParameter<double>("HTMin")), NjMin(iConfig.getParameter<int>("NjMin")),
  EleVetoIdMapToken(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap2016"))),
  EleLooseIdMapToken(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap2016"))),
  EleMediumIdMapToken(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap2016"))),
  EleTightIdMapToken(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap2016")))
  //mvaValuesMapToken(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
  //mvaCategoriesMapToken(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))
{
   //Explicitly declare shared resource TFileService to make it threadsafe
   usesResource("TFileService");
   counter = 0;

   ////////////////////
   ////// Tokens //////
   ////////////////////
   JetToken = consumes<std::vector<pat::Jet> >(edm::InputTag("selectedUpdatedPatJetsDeepCSV"));
   MuonToken = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
   ElectronToken = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons"));
   GsfElectronToken = consumes<edm::View<reco::GsfElectron> >(edm::InputTag("slimmedElectrons"));
   METToken = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));
   VtxToken = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));
   HLTToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
   FltToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO"));
   //ConversionsToken = consumes<std::vector<reco::Conversion> >(edm::InputTag("reducedConversions")); //EleID has conversion veto
   BeamSpotToken = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
   //RhoToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll")); //not necessary at this moment
   GenToken = consumes<std::vector<pat::PackedGenParticle> >(edm::InputTag("packedGenParticles"));
   PrunedGenToken = consumes<std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
   GenJetToken = consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets"));
   GenJetAK8Token = consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJetsAK8"));

   //MET Filter Settings
   //HBHENoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("HBHENoiseFilter_Selector_");
   //EEBadScNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("EEBadScNoiseFilter_Selector_");
   //HBHENoiseFilter_Selector_ = "HBHENoiseFilter_Selector_";
   //EEBadScNoiseFilter_Selector_ = "EEBadScNoiseFilter_Selector_";

   ///////////////////////////////
   /// HT Minimum Notification ///
   ///////////////////////////////
   std::cout << "The Minimum HT Cut is currently: " << HTMin << "GeV. Set HTMin = <value> in python configuration of ntupler." << std::endl;
   //////////////////////////////////////////////
   /// Per-Year definitions for future ReReco ///
   //////////////////////////////////////////////

   if(is2016){
   ///////////
   /// HLT ///
   ///////////
   //Store the string name of triggers in arrays with postfix _S, and corresponding booleans will be in arrays with postfix _B, int representation of bits with postifx _Bits
     std::cout << "Defining triggers for 2016" << std::endl;
     //MuMu Triggers
     HLT_MuMu_S.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
     HLT_MuMu_S.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
     HLT_MuMu_S.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
     HLT_MuMu_S.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

     //ElMu Triggers
     HLT_ElMu_S.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
     HLT_ElMu_S.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
     HLT_ElMu_S.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
     HLT_ElMu_S.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");

     //ElEl Triggers
     HLT_ElEl_S.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

     //Mu Triggers
     HLT_Mu_S.push_back("HLT_IsoTkMu24_v");
     HLT_Mu_S.push_back("HLT_IsoMu24_v");

     //El Triggers
     HLT_El_S.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");

     //MET Filters
     MET_Flt_S.push_back("FIXME"); //not used for MC, events failing filters aren't kept

   }
   else if(is2017){
     std::cout << "Defining triggers for 2017" << std::endl;
     std::cout << "FIXME!" << std::endl;
     throw cms::Exception("2017 triggers not defined!");
   }
   else if(is2018){
     std::cout << "Defining triggers for 2018" << std::endl;
     std::cout << "FIXME!" << std::endl;
     throw cms::Exception("2018 triggers not defined!");
   }
   else{
     std::cout << "Error: Data is not from 2016, 2017, or 2018. Is it ReReco?" << std::endl;
     throw cms::Exception("Undefined data year (2016-2018)!");
   }

   //Set HLT bits to zero via initialization, using the boost::dynamic_bitset, with the proper number of bits
   HLT_MuMu_B = boost::dynamic_bitset<>(HLT_MuMu_S.size(), 0ul); //sets number of bits corresponding to triggers per channel, and 0 unsigned long value
   HLT_ElMu_B = boost::dynamic_bitset<>(HLT_ElMu_S.size(), 0ul);
   HLT_ElEl_B = boost::dynamic_bitset<>(HLT_ElEl_S.size(), 0ul);
   HLT_Mu_B = boost::dynamic_bitset<>(HLT_Mu_S.size(), 0ul);
   HLT_El_B = boost::dynamic_bitset<>(HLT_El_S.size(), 0ul);

   //Set MET bits to zero via initialization
   MET_Flt_B = boost::dynamic_bitset<>(MET_Flt_S.size(), 0ul);

}


SLntupler::~SLntupler()
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
SLntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   counter++;
   if(deBug) std::cout << "Counter: " << counter << std::endl;

   using namespace edm;
   edm::Handle<std::vector<pat::Jet> > jets;
   iEvent.getByToken(JetToken, jets);
   if(!jets.isValid()) {
     throw cms::Exception("Jet collection not valid!"); 
   }
   edm::Handle<std::vector<pat::Muon> > muons;
   iEvent.getByToken(MuonToken, muons);
   if(!muons.isValid()) {
     throw cms::Exception("Muon collection not valid!"); 
   }
   // edm::Handle<std::vector<pat::Electron> > electrons;
   // iEvent.getByToken(ElectronToken, electrons);
   // if(!electrons.isValid()) {
   //   throw cms::Exception("Electron collection not valid!"); 
   // }
   edm::Handle<edm::View<reco::GsfElectron> > gsfelectrons; //Just a cast of pat::Electrons, so we can do gymnastics for IDs... not necessary in later CMSSW releases
   iEvent.getByToken(GsfElectronToken, gsfelectrons);
   if(!gsfelectrons.isValid()) {
     throw cms::Exception("Casting electrons to reco::GsfElectron failed!");
   }
   edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
   iEvent.getByToken(EleVetoIdMapToken,veto_id_decisions);
   if(!veto_id_decisions.isValid()){
     throw::cms::Exception("Ele Cut-based Veto ID decisions not valid!");    // Note:  VID ID modules must have been run upstream.
   }
   edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
   iEvent.getByToken(EleLooseIdMapToken,loose_id_decisions);
   if(!loose_id_decisions.isValid()){
     throw::cms::Exception("Ele Cut-based Loose ID decisions not valid!");    // Note:  VID ID modules must have been run upstream.
   }
   edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   iEvent.getByToken(EleMediumIdMapToken,medium_id_decisions);
   if(!medium_id_decisions.isValid()){
     throw::cms::Exception("Ele Cut-based Medium ID decisions not valid!");    // Note:  VID ID modules must have been run upstream.
   }
   edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
   iEvent.getByToken(EleTightIdMapToken,tight_id_decisions);
   if(!tight_id_decisions.isValid()){
     throw::cms::Exception("Ele Cut-based Tight ID decisions(s) not valid!");    // Note:  VID ID modules must have been run upstream.
   }
   // Get MVA values and categories (optional)
   // edm::Handle<edm::ValueMap<float> > mvaValues;
   // edm::Handle<edm::ValueMap<int> > mvaCategories;
   // iEvent.getByToken(mvaValuesMapToken_,mvaValues);
   // iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);
   edm::Handle<std::vector<pat::MET> > mets;
   iEvent.getByToken(METToken, mets);
   if(!mets.isValid()) {
     throw cms::Exception("MET collection not valid!"); 
   }
   edm::Handle<edm::TriggerResults> METFlt;
   iEvent.getByToken(FltToken, METFlt);
   if(!METFlt.isValid()) {
     throw cms::Exception("MET Filter collection not valid!"); 
   }
   edm::Handle<edm::TriggerResults> HLTTrg;
   iEvent.getByToken(HLTToken, HLTTrg);
   if(!HLTTrg.isValid()) {
     throw cms::Exception("HLT collection not valid!"); 
   }
   edm::Handle<std::vector<reco::Vertex> > vertices;
   iEvent.getByToken(VtxToken, vertices);
   if(!vertices.isValid()) {
     throw cms::Exception("Vertex collection not valid!"); 
   } 
   // edm::Handle<reco::ConversionCollection> conversions;
   // iEvent.getByToken(ConversionsToken, conversions);
   // if(!conversions.isValid()) {
   //   throw cms::Exception("Conversions collection not valid!");
   // }
   edm::Handle<reco::BeamSpot> theBeamSpot;
   iEvent.getByToken(BeamSpotToken,theBeamSpot); 
   if(!theBeamSpot.isValid()) {
     throw cms::Exception("BeamSpot not valid!");
   }
   // edm::Handle<double> rho;
   // iEvent.getByToken(RhoToken, rho);
   //protect against bad parameter?
   if(isMC){
     edm::Handle <std::vector<pat::PackedGenParticle> > gens;
     iEvent.getByToken(GenToken, gens);
     if(!gens.isValid()) {
       throw cms::Exception("Gen Particle collection not valid!");
     }

     edm::Handle <std::vector<reco::GenParticle> > prunedGens;
     iEvent.getByToken(PrunedGenToken, prunedGens);
     if(!prunedGens.isValid()) {
       throw cms::Exception("Pruned Gen Particle collection not valid!");
     }

     edm::Handle <std::vector<reco::GenJet> > genjets;
     iEvent.getByToken(GenJetToken, genjets);
     if(!prunedGens.isValid()) {
       throw cms::Exception("Gen Jet collection not valid!");
     }

     edm::Handle <std::vector<reco::GenJet> > genjetsak8;
     iEvent.getByToken(GenJetAK8Token, genjetsak8);
     if(!prunedGens.isValid()) {
       throw cms::Exception("Gen Jet AK8 collection not valid!");
     }
   }

   //See /afs/cern.ch/user/n/nmangane/DAS/EGammaExercises/MuonExercise3/plugins/MuonExercise3.cc for histogram array creation, tight (vertex) selection
   //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ CP code
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


   //Clear pointers
   JetLVec->clear();
   selectedLepLVec->clear();
   METLVec->clear();
   VertexVec->clear();
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

   //reset bitsets
   HLT_MuMu_B.reset();
   HLT_ElMu_B.reset();
   HLT_ElEl_B.reset();
   HLT_Mu_B.reset();
   HLT_El_B.reset();
   MET_Flt_B.reset();

   /////////////////////////////
   /// HLT TRIGGER SELECTION ///
   /////////////////////////////
   if(verBose) std::cout << "\n=============HLT=============\n";
   const edm::TriggerNames &HLTnames = iEvent.triggerNames(*HLTTrg);
   for (unsigned int i = 0, n = HLTTrg->size(); i < n; ++i) {
     //std::cout << HLTnames.triggerName(i) << std::endl;
     //FIXME Option: Create a mapping on a per-year basis to try matching based on position, and if that fails, then loop through instead
     //For efficicency: Do one search through names to find locations, then store them in global variables to be used in all following events
     std::string trgName = HLTnames.triggerName(i); //string of triggername
     bool trgBit = HLTTrg->accept(i); //accept bit
     std::string delimeter = "_v"; //common delimeter in triggers
     std::string trgSubName = trgName.substr(0, trgName.find(delimeter) + delimeter.length()); //only take the common portion of trigger (up through version marker _v)
     for(uint j = 0; j < HLT_MuMu_S.size(); j++)
       if (trgSubName == HLT_MuMu_S[j]){
	 if(verBose) std::cout << " Name: " << trgName << " Initial Bits: " << HLT_MuMu_B;
	 HLT_MuMu_B[j] = trgBit; //sets individual bit, starting from most significant ("leftmost" in 'operator<<' language)
	 if(verBose) std::cout << " Accepted: " << trgBit << " Final Bits: " << HLT_MuMu_B << std::endl;
       }
     for(uint jj = 0; jj < HLT_ElMu_S.size(); jj++)
       if (trgSubName == HLT_ElMu_S[jj]){
	 if(verBose) std::cout << " Name: " << trgName << " Initial Bits: " << HLT_ElMu_B;
	 HLT_ElMu_B[jj] = trgBit;
	 if(verBose) std::cout << " Accepted: " << trgBit << " Final Bits: " << HLT_ElMu_B << std::endl;
       }
     for(uint jjj = 0; jjj < HLT_ElEl_S.size(); jjj++)
       if (trgSubName == HLT_ElEl_S[jjj]){
	 if(verBose) std::cout << " Name: " << trgName << " Initial Bits: " << HLT_ElEl_B;
	 HLT_ElEl_B[jjj] = trgBit;
	 if(verBose) std::cout << " Accepted: " << trgBit << " Final Bits: " << HLT_ElEl_B << std::endl;
       }
     for(uint k = 0; k < HLT_Mu_S.size(); k++)
       if (trgSubName == HLT_Mu_S[k]){
	 if(verBose) std::cout << " Name: " << trgName << " Initial Bits: " << HLT_Mu_B;
	 HLT_Mu_B[k] = trgBit;
	 if(verBose) std::cout << " Accepted: " << trgBit << " Final Bits: " << HLT_Mu_B << std::endl;
       }
     for(uint kk = 0; kk < HLT_El_S.size(); kk++)
       if (trgSubName == HLT_El_S[kk]){
	 if(verBose) std::cout << " Name: " << trgName << " Initial Bits: " << HLT_El_B;
	 HLT_El_B[kk] = trgBit;
	 if(verBose) std::cout << " Accepted: " << trgBit << " Final Bits: " << HLT_El_B << std::endl;
       }     
   }
   //bitset -> unsigned long -> unsigned int (ROOT-compatible format)
   HLT_MuMu_Bits = HLT_MuMu_B.to_ulong(); 
   HLT_ElMu_Bits = HLT_ElMu_B.to_ulong();
   HLT_ElEl_Bits = HLT_ElEl_B.to_ulong();
   HLT_Mu_Bits = HLT_Mu_B.to_ulong();
   HLT_El_Bits = HLT_El_B.to_ulong();

   //SL + DL Triggers
   // if(!(HLT_MuMu_B.any() || HLT_ElMu_B.any() || HLT_ElEl_B.any() || HLT_Mu_B.any() || HLT_El_B.any() ) ) //if NO triggers pass, want to skip rest of this module
   //   return; //void main, so just return nothing 

   //SL Trigger only
   if(!(HLT_Mu_B.any() || HLT_El_B.any())){
     if(verBose)
       std::cout << "No SL triggers pass!" << std::endl;
     return; //move to next event if not SL trigger
   }
   //===/////////////
   //===//// CUT ///
   //===///////////

   ////////////////////////////
   /// MET FILTER SELECTION ///
   ////////////////////////////
   if(isData){
     const edm::TriggerNames &names = iEvent.triggerNames(*METFlt);
     for (uint i = 0; i < METFlt->size(); ++i) {
       std::string fltName = names.triggerName(i); //convenient name storage
       if(verBose) std::cout << fltName << std::endl;
       bool fltBit = METFlt->accept(i); //filter pass bit
       for(uint j = 0; j < MET_Flt_S.size(); j++){
	 if (fltName == MET_Flt_S[j]){
	   if(verBose) std::cout << "Initial Bits: " << MET_Flt_B;
	   MET_Flt_B[j] = fltBit; //store bit decision in bitset
	   if(verBose) std::cout << " Name: " << fltName << " Accepted: " << fltBit << " Bits: " << MET_Flt_B << std::endl;
	 }     
       }
     }
   }
   MET_Flt_Bits = MET_Flt_B.to_ulong(); //store flt decision even if default 0
   // End filters stuff


   ////////////////////
   //// Event info ////
   ////////////////////
   nRun = iEvent.id().run();
   nLumiBlock = iEvent.id().luminosityBlock();
   nEvent = iEvent.id().event();

   /////////////////////
   //// Good Vertex ////
   /////////////////////
   std::vector<reco::Vertex>::const_iterator firstGoodVertex = vertices->end();

   for (std::vector<reco::Vertex>::const_iterator verts=vertices->begin(); verts!=firstGoodVertex; verts++) {
     if (!verts->isFake() && verts->ndof()>4 && verts->position().Rho()<2. && std::abs(verts->position().Z())<24.) {
       if(firstGoodVertex == vertices->end()) firstGoodVertex = verts;
       break;
     }
   }
   //Require good vertex
   if(firstGoodVertex == vertices->end()){
     if(verBose)
       std::cout << "No good PV!" << std::endl;
     return;
   }
   TLorentzVector perVertexLVec;
   perVertexLVec.SetXYZT(firstGoodVertex->position().X(), firstGoodVertex->position().Y(), firstGoodVertex->position().Z(), 0); //dummy 0 for Time coordinate
   VertexVec->push_back(perVertexLVec); //First vertex in vector is primary. SVs can be emplaced afterwards
   //===/////////////
   //===//// CUT ///
   //===///////////

   ///////////////////////
   //// MET Selection ////
   ///////////////////////
   if(verBose) std::cout << "\n=============MET=============\n";
   const pat::MET &met = mets->front();
   if(verBose) std::cout << "MET Pt: " << met.pt() << " Phi: " << met.phi() << std::endl;
   if(met.pt() < 50){
     if(verBose) 
       std::cout << "MET below 50GeV!" << std::endl;
     return;
   }
   // TLorentzVector perMETLVec;
   // perMETLVec.SetPtEtaPhiE(met.pt(), met.eta(), met.phi(), met.energy());
   // METLVec->push_back(perMETLVec);
   //===/////////////
   //===//// CUT ///
   //===///////////

   ////////////////////////
   //// Selected Muons ////
   ////////////////////////
   if(verBose) std::cout << "\n============Muons============\n";
   for(const pat::Muon& muon : *muons){
     //Min cuts for the Loose (veto) Muon selection
     if(muon.pt() <= 10 || fabs(muon.eta()) >= 2.5)
       continue;
   //===/////////////
   //===//// CUT ///
   //===///////////

     //Calculate the Relative Isolation
     double relIso = (muon.pfIsolationR04().sumChargedHadronPt + fmax(0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt))/muon.pt();

     //Cut on max isolation for the loose(veto) muon
     if(relIso >= 0.25)
       continue;
   //===/////////////
   //===//// CUT ///
   //===///////////

     //Select Tight Muons (Tight ID, relIso < 0.15, pt > 25, |eta| < 2.1)
     if(muon.isTightMuon(*firstGoodVertex) && muon.pt() > 25 && fabs(muon.eta()) < 2.1 && relIso < 0.15 ){
       //setup LVec for Muon
       TLorentzVector perMuonLVec; 
       perMuonLVec.SetPtEtaPhiE( muon.pt(), muon.eta(), muon.phi(), muon.energy() );

       //add to selected leptons
       selectedLepLVec->push_back(perMuonLVec);

       //debug info
       if(verBose)
	 std::cout << "Pt: " << muon.pt() << " Eta: " << muon.eta() << " Phi: " << muon.phi()  << "Tight ID: " << muon.isTightMuon(*firstGoodVertex)
		   << " Loose ID: " << muon.isLooseMuon() << " relIso: " << relIso << std::endl;
     }

     //Select Loose Muons for Veto (Loose ID, relIso < 0.25, pt > 10, |eta| < 2.5)
     else if(muon.isLooseMuon() ){
       if(verBose)
	 std::cout << "Loose Muon (2nd Lepton) detected!" << std::endl;
       return;  // If there are any Loose/Veto Leptons, then not in the SL channel!
     }
   //===/////////////
   //===//// CUT ///
   //===///////////
   }
   

   ////////////////////////////
   //// Selected Electrons ////
   ////////////////////////////
   if(verBose) std::cout << "\n==========Electrons==========\n";
   int ii = -1;
   for(const reco::GsfElectron& electron : *gsfelectrons){
     ii++;
     //Select veto lepton cuts minimum
     if(electron.pt() < 15 || fabs(electron.eta()) > 2.5 )
       continue;
   //===/////////////
   //===//// CUT ///
   //===///////////

     //Need to do fucking juggling to get a GsfElectron pointer so that we can then get the ID. Thankfully fixed in MiniAOD V2 (94X +), where
     //electron ID is just accessed with electronIterator->electronID("Name_Of_The_ID")
     const auto el = gsfelectrons->ptrAt(ii);
     bool vetoID = (*veto_id_decisions)[el];
     bool tightID = (*tight_id_decisions)[el];

     //verBose info
     if(verBose)
       std::cout << "Pt: " << electron.pt() << " Eta: " << electron.eta() << " Phi: " << electron.phi() << " Veto ID: " << vetoID << " Tight ID: " << tightID << std::endl;

     //Calculate Isolation
     //https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/plugins/IsoValueMapProducer.cc
     // auto iso = electron.pfIsolationVariables();
     // auto chg = iso.sumChargedHadronPt;
     // auto neu = iso.sumNeutralHadronEt;
     // auto pho = iso.sumPhotonEt;
     //auto ea = ea_pfiso_->getEffectiveArea(fabs(getEtaForEA(&electron)));
     //float scale = relative_ ? 1.0/obj.pt() : 1;

     if(electron.pt() > 30 && fabs(electron.eta()) < 2.1 && tightID){
     //Set up Lorentz Vector for Electrons
     TLorentzVector perElectronLVec;
     perElectronLVec.SetPtEtaPhiE( electron.pt(), electron.eta(), electron.phi(), electron.energy() );

     //Add to selected leptons
     selectedLepLVec->push_back(perElectronLVec);
     }
     else if(vetoID){
       if(verBose)
	 std::cout << "Veto electron (2nd Lepton) detected!" << std::endl;
       return;
     }
   //===/////////////
   //===//// CUT ///
   //===///////////
   }
       
   //SL Select events with only 1 isolated lepton
   if(selectedLepLVec->size() > 1){
     if(verBose)
       std::cout << "Two well-ID'd, well-isolated leptons detected!" << std::endl;
     return;
   }
   if(selectedLepLVec->size() < 1){
     if(verBose)
       std::cout << "No well-ID'd, well isolated leptons detected!" << std::endl;
     return;
   }
   //===/////////////
   //===//// CUT ///
   //===///////////

   ///////////////////////
   //// Selected Jets ////
   ///////////////////////
   if(verBose) std::cout << "\n=============Jets============\n";
   //Reset HT to 0
   HT = 0;
   if(deBug && counter == theProblemEvent) std::cout << "L1 ";
   for(const pat::Jet&jet : *jets){
     //Jet Selection 30GeV, usual calorimeter acceptance
     if(jet.pt() < 30 || fabs(jet.eta()) >= 2.5)
       continue;
   //===/////////////
   //===//// CUT ///
   //===///////////
   if(deBug && counter == theProblemEvent) std::cout << "L2 ";     
     TLorentzVector perJetLVec;
     perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
     double qgPtD = jet.userFloat("QGTagger:ptD");
     double qgAxis1 = jet.userFloat("QGTagger:axis1");
     double qgAxis2 = jet.userFloat("QGTagger:axis2");
     double qgMult = static_cast<double>(jet.userInt("QGTagger:mult"));
     double deepCSVb;
     double deepCSVc;
     double deepCSVl;
     double deepCSVbb;
     double deepCSVcc;
     if(maskDeepCSV){
       deepCSVb = -1000;
       deepCSVc = -1000; 
       deepCSVl = -1000;
       deepCSVbb = -1000;
       deepCSVcc = -1000;
     }
     else{
       deepCSVb = jet.bDiscriminator("pfDeepCSVJetTags:probb");
       deepCSVc = jet.bDiscriminator("pfDeepCSVJetTags:probc");
       deepCSVl = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
       deepCSVbb = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
       deepCSVcc = jet.bDiscriminator("pfDeepCSVJetTags:probcc");
     }
     if(deBug && counter == theProblemEvent) std::cout << "L3 ";
     double btag = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
     double chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
     double neutralHadronEnergyFraction = jet.neutralHadronEnergyFraction();
     double chargedEmEnergyFraction = jet.chargedEmEnergyFraction();
     double neutralEmEnergyFraction = jet.neutralEmEnergyFraction();
     double muonEnergyFraction = jet.muonEnergyFraction();
     double photonEnergyFraction = jet.photonEnergyFraction();
     double electronEnergyFraction = jet.electronEnergyFraction();
   if(deBug && counter == theProblemEvent) std::cout << "L4 ";
     double recoJetsHFHadronEnergyFraction = jet.HFHadronEnergyFraction();
     double recoJetsHFEMEnergyFraction = jet.HFEMEnergyFraction();
     double chargedHadronMultiplicity = jet.chargedHadronMultiplicity();
     double neutralHadronMultiplicity = jet.neutralHadronMultiplicity();
     double photonMultiplicity = jet.photonMultiplicity();
     double electronMultiplicity = jet.electronMultiplicity();
     double muonMultiplicity = jet.muonMultiplicity();
   if(deBug && counter == theProblemEvent) std::cout << "L5 ";
     //Jet ID Loose selection for 2016 (doesn't exist for 2017 or 2018! Tight selection must be used!)
     bool looseJetID = 
       (neutralHadronEnergyFraction < 0.99 && neutralEmEnergyFraction < 0.99 && (jet.chargedMultiplicity() + jet.neutralMultiplicity() ) > 1) && 
       ((fabs(jet.eta()) <= 2.4 && chargedHadronEnergyFraction > 0 && jet.chargedMultiplicity() > 0 && chargedEmEnergyFraction < 0.99) || fabs(jet.eta()) > 2.4 );
     if(verBose) std::cout << "Pt: " << jet.pt() << " Eta: " << jet.eta() << " Phi: " << jet.phi() << " Jet Loose ID: " 
			 << looseJetID << " CSVv2: " << btag << " DeepCSV(b+bb): " << deepCSVb + deepCSVbb << std::endl;
   if(deBug && counter == theProblemEvent) std::cout << "L6 ";
     if(!looseJetID)
       continue;
   //===/////////////
   //===//// CUT ///
   //===///////////
   if(deBug && counter == theProblemEvent) std::cout << "L7 ";
   // try{
   //   selectedLepLVec->at(0);
   // }
   // catch{
   //   std::cout << "It appears, good Madam/Sir, that this Lepton Vector is indeed empty! Hats off!" << std::endl;
   // }
     double dR = perJetLVec.DeltaR(selectedLepLVec->at(0));
   if(deBug && counter == theProblemEvent) std::cout << "L7B ";
     if(verBose) 
       std::cout << ">>Cross-Cleaning<< Jet Eta: " << jet.eta() << " Jet Phi: " << jet.phi() << " Lep Eta: " << selectedLepLVec->at(0).Eta() 
		 << " Lep Phi: " << selectedLepLVec->at(0).Phi() << " Jet-Lep DeltaR: " << dR << std::endl;
     //0.4 for SL, 0.3 for DL!
     if(dR < 0.4){
       if(verBose)
	 std::cout << "Cross-cleaning jet!" << std::endl;
       continue;
     }
   if(deBug && counter == theProblemEvent) std::cout << "L8 ";
   //===/////////////
   //===//// CUT ///
   //===///////////

     //Sum HT for selected jets
     HT += jet.pt();

   if(deBug && counter == theProblemEvent) std::cout << "L9 ";
     JetLVec->push_back(perJetLVec);
     qgPtDVec->push_back(qgPtD);
     qgAxis1Vec->push_back(qgAxis1);
     qgAxis2Vec->push_back(qgAxis2);
     qgMultVec->push_back(qgMult);
   if(deBug && counter == theProblemEvent) std::cout << "L10 ";
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
   if(deBug && counter == theProblemEvent) std::cout << "L11 ";
     electronEnergyFractionVec->push_back(electronEnergyFraction);
     recoJetsHFHadronEnergyFractionVec->push_back(recoJetsHFHadronEnergyFraction);
     recoJetsHFEMEnergyFractionVec->push_back(recoJetsHFEMEnergyFraction);
     chargedHadronMultiplicityVec->push_back(chargedHadronMultiplicity);
     neutralHadronMultiplicityVec->push_back(neutralHadronMultiplicity);
     photonMultiplicityVec->push_back(photonMultiplicity);
     electronMultiplicityVec->push_back(electronMultiplicity);
     muonMultiplicityVec->push_back(muonMultiplicity);
   if(deBug && counter == theProblemEvent) std::cout << "L12 ";
   }

   //////////////////////
   /// HT Minimum cut ///
   //////////////////////
   if(verBose)
     std::cout << "HT: " << HT << std::endl;
   if(HT < HTMin){
     if(verBose)
       std::cout << "HT below threshold of " << HTMin << "!" << std::endl;
     return;
   }
   if(deBug && counter == theProblemEvent) std::cout << "L13 ";
   //===/////////////
   //===//// CUT ///
   //===///////////

   ////////////////////////
   /// Njet Minimum cut ///
   ////////////////////////
   if(verBose)
     std::cout << "Njet: " << JetLVec->size() << std::endl;
   if(NjMin > -1)
     if(JetLVec->size() < (uint)NjMin){
       if(verBose) 
	 std::cout << "Below threshold of minimum jets, continuing to next event" << std::endl;
       return;
     }
   if(deBug && counter == theProblemEvent) std::cout << "L14 ";
   nEvts++;

   ///////////////////////////////////////////
   /// Fill Tree (written by TFileService) ///
   ///////////////////////////////////////////
   tree->Fill();
   if(deBug && counter == theProblemEvent) std::cout << "L15 ";

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
SLntupler::beginJob()
{   
   nEvts = 0;
   nRun = -1;
   nLumiBlock = -1; 
   nEvent = -1;
   HT = -1;
   //FIXME: Add missing variables for leptons, isolation, jetID, HT, etc.
   JetLVec = new std::vector<TLorentzVector>;
   METLVec = new std::vector<TLorentzVector>;
   selectedLepLVec = new std::vector<TLorentzVector>;
   vetoLepLVec = new std::vector<TLorentzVector>;
   VertexVec = new std::vector<TLorentzVector>;
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
   tree = fs->make<TTree>("nTuple", "Event nTuple");
   tree->Branch("run", &nRun);
   tree->Branch("lumiBlock", &nLumiBlock);
   tree->Branch("event", &nEvent);
   tree->Branch("HT", &HT);
   // tree->Branch("HTX", &HTX);
   // tree->Branch("HT2M", &HT2M);
   tree->Branch("HLT_MuMu_Bits", &HLT_MuMu_Bits);
   tree->Branch("HLT_ElMu_Bits", &HLT_ElMu_Bits);
   tree->Branch("HLT_ElEl_Bits", &HLT_ElEl_Bits);
   tree->Branch("HLT_Mu_Bits", &HLT_Mu_Bits);
   tree->Branch("HLT_El_Bits", &HLT_El_Bits);
   tree->Branch("MET_Flt_Bits", &MET_Flt_Bits);
   tree->Branch("JetLVec", "vector<TLorentzVector>", &JetLVec, 32000,-1);
   tree->Branch("METLVec", "vector<TLorentzVector>", &METLVec, 32000,-1);
   tree->Branch("selectedLepLVec", "vector<TLorentzVector>", &selectedLepLVec, 32000,-1);
   tree->Branch("VertexVec", "vector<TLorentzVector>", &VertexVec, 32000, -1);
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
SLntupler::endJob() 
{
   usesResource("TFileService");
   //********** WARNING **********//
   //Using tree->Write(***) will invoke Write() TWICE, as the cmsRun framework does this for you. This results in a duplicate tree in the file
   // tree->GetDirectory()->cd();
   // tree->Write("", TObject::kOverwrite);

   //Clear pointers
   JetLVec->clear();
   selectedLepLVec->clear();
   METLVec->clear();
   vetoLepLVec->clear();
   VertexVec->clear();
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

   delete JetLVec;
   delete selectedLepLVec;
   delete METLVec;
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
SLntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SLntupler);
