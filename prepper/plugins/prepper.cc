// -*- C++ -*-
//
// Package:    Demo/prepper
// Class:      prepper
// 
/**\class prepper prepper.cc Demo/prepper/plugins/prepper.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nick Manganelli
//         Created:  Thu, 23 Aug 2018 21:15:14 GMT
//
//


// system include files
//#include <stdio>
//#include "TTree.h"
//#include "TNtuple.h"
//#include <TMatrixDSym.h>
//#include <TMatrixDSymEigen.h>
//#include <TVectorD.h>
//#include <ctime>
//#include <fstream>
//#include <sstream>

#include <memory>
#include <cmath>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

//#include "CommonTools/UtilAlgos/interface/TFileService.h"


// #include "FWCore/Framework/interface/Event.h"
// #include "FWCore/Framework/interface/MakerMacros.h"

// #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/Utilities/interface/StreamID.h"

// #include "DataFormats/Common/interface/Handle.h"
// #include "DataFormats/Common/interface/Ref.h"
// #include "DataFormats/Common/interface/RefToBase.h"
// #include "DataFormats/Common/interface/RefVector.h"
// #include "DataFormats/Common/interface/ValueMap.h"
// #include "DataFormats/Candidate/interface/CandidateFwd.h"
// #include "DataFormats/Candidate/interface/CandMatchMap.h"
// #include "DataFormats/JetReco/interface/PFJet.h"
// #include "DataFormats/JetReco/interface/GenJet.h"
// #include "DataFormats/JetReco/interface/JetCollection.h"
// #include "DataFormats/PatCandidates/interface/Jet.h"
// #include "DataFormats/PatCandidates/interface/Tau.h"
// #include "DataFormats/MuonReco/interface/MuonFwd.h"
// #include "DataFormats/PatCandidates/interface/Muon.h"
// #include "DataFormats/PatCandidates/interface/Electron.h"
// #include "DataFormats/PatCandidates/interface/MET.h"
// #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
// #include "DataFormats/Math/interface/deltaR.h"
// #include "DataFormats/Math/interface/deltaPhi.h"
// #include "DataFormats/VertexReco/interface/Vertex.h"
// #include "DataFormats/VertexReco/interface/VertexFwd.h"
// #include "DataFormats/TrackReco/interface/TrackFwd.h"


//
// class declaration
//

class prepper : public edm::stream::EDProducer<> {
   public:
      explicit prepper(const edm::ParameterSet&);
      ~prepper();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<pat::Jet>> JetToken;
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
prepper::prepper(const edm::ParameterSet& iConfig)
{
   JetToken = consumes<std::vector<pat::Jet> >(edm::InputTag("selectedUpdatedPatJetsDeepCSV"));
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


prepper::~prepper()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
prepper::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
   //std::cout << "Creating Handle..." << std::endl;
   edm::Handle<std::vector<pat::Jet> > jets;
   //std::cout << "grabbing product by token..." << std::endl;
   iEvent.getByToken(JetToken, jets);
   //std::cout << "Collection acquired" << std::endl;
   int i = 0;
   for(const pat::Jet& jet : *jets)
   //   for(std::vector<pat::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet)
     {
       i++;
       //No need to keep jets below 20 GeV
       if(jet.pt() < 20) continue;

       TLorentzVector perJetLVec;
       perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );

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
       
       std::cout << " PtD: " << qgPtD << " Ax1: " << qgAxis1 << " Ax2: " << qgAxis2 << " QG_Mult: " << qgMult << std::endl;
       std::cout << " Pt: " << jet.pt() << " DeepCSVb+bb: " << (deepCSVb + deepCSVbb) << " HFHadEnFrac: " << recoJetsHFHadronEnergyFraction << std::endl;
       }
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
prepper::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
prepper::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
prepper::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
prepper::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
prepper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
prepper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
prepper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(prepper);
