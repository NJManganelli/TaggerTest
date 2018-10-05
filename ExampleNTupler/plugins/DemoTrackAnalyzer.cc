// -*- C++ -*-
//
// Package:    Demo/DemoTrackAnalyzer
// Class:      DemoTrackAnalyzer
// 
/**\class DemoTrackAnalyzer DemoTrackAnalyzer.cc Demo/DemoTrackAnalyzer/plugins/DemoTrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jhovanny Andres Mejia Guisao
//         Created:  Tue, 05 Sep 2017 22:37:14 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"



//
// class declaration
//
class DemoTrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoTrackAnalyzer(const edm::ParameterSet&);
      ~DemoTrackAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackTags_;
      TH1D *histo;
      TH1D *histo_dxy;
      TH1D *histo_dxyError;

      TH1D *histo_Hits;
      TH1D *histo_PixelHits;

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
DemoTrackAnalyzer::DemoTrackAnalyzer(const edm::ParameterSet& iConfig) :
trackTags_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("tracks")))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  histo = fs->make<TH1D>("pt" , "Pt" , 50 , 0 , 50 );

  //other histos
  histo_dxy = fs->make<TH1D>("dxy" , "dxy" , 50 , 0 , 10 );
  histo_dxyError = fs->make<TH1D>("dxyError" , "dxyError" , 50 , 0 , 1 );

  histo_Hits = fs->make<TH1D>("Hits" , "Hits" , 30 , 0 , 30 );
  histo_PixelHits = fs->make<TH1D>("PixelHits" , "PixelHits" , 10 , 0 , 10 );
}


DemoTrackAnalyzer::~DemoTrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace std;
    
  

   //Handle<TrackCollection> tracks;
   //iEvent.getByLabel(trackTags_,tracks);
   edm::Handle< View<pat::PackedCandidate> > tracks;
   iEvent.getByToken(trackTags_,tracks);

   //std::cout<< " number of tracks:   " << tracks->size() << endl;
   //LogInfo("Demo") << "number of tracks "<<tracks->size();

   for(View<pat::PackedCandidate>::const_iterator itTrack = tracks->begin();
       itTrack != tracks->end(); ++itTrack )
     {

       //if you want to use dxy or dz you need to be sure the pt
       //of the tracks is bigger than 0.5 GeV, otherwise you will get error

       if(itTrack->pt()<=0.5)continue; 
       if(itTrack->charge()==0) continue;// NO neutral objects
       if(fabs(itTrack->pdgId())!=211) continue;//Due to the lack of the particle ID all the tracks for cms are pions(ID==211)
       if(!(itTrack->trackHighPurity())) continue;

       
       //std::cout<< " track momentum components:   " << itTrack->pt()<<", "<<itTrack->px()<<", "<<itTrack->py()<<", "<<itTrack->pz() << endl;
       //std::cout<< " track impact parameter: " << itTrack->dz()<<", "<<itTrack->dzError()<<itTrack->dxy()<<", "<<itTrack->dxyError()< endl;
       //std::cout<< " others track properties: " <<itTrack->eta()<<", "<<itTrack->phi() <<", "<<itTrack->charge() << endl;

       
       // **** fill Pt histogram
       histo->Fill(itTrack->pt());

      
       //*****Accessing track impact parameters **********
       histo_dxy->Fill(itTrack->dxy());
       histo_dxyError->Fill(itTrack->dxyError());


       //****** Accessing the hit pattern of a track  ********
       /*
       // NOTA BENE: please be carefully. From MiniAOD link you will se the next comment:
       " approximate hitPattern() and trackerExpectedHitsInner() that yield the correct number of hits, pixel hits and the information returned by los        tInnerHits()
       (but nothing more, so e.g. you can't rely on methods counting layers, or expected outer hits, or lost hits within the track) "
       */
             
       histo_Hits->Fill(itTrack->numberOfHits());
       histo_PixelHits->Fill(itTrack->numberOfPixelHits());
       

       
     }

   /*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   */
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
DemoTrackAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoTrackAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
DemoTrackAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoTrackAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoTrackAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoTrackAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoTrackAnalyzer);
