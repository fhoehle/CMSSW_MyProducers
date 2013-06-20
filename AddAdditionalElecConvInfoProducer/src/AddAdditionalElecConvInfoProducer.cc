// -*- C++ -*-
//
// Package:    AddAdditionalElecConvInfoProducer
// Class:      AddAdditionalElecConvInfoProducer
// 
/**\class AddAdditionalElecConvInfoProducer AddAdditionalElecConvInfoProducer.cc CMSSW_MyProducers/AddAdditionalElecConvInfoProducer/src/AddAdditionalElecConvInfoProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Felix Hoehle
//         Created:  Thu Oct  4 16:58:14 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
////#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
////#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
////#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
////#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
////#include "DataFormats/TrackReco/interface/Track.h"
////#include "DataFormats/TrackReco/interface/TrackFwd.h"
////#include "DataFormats/Scalers/interface/DcsStatus.h"
//
// class declaration
//

class AddAdditionalElecConvInfoProducer : public edm::EDProducer {
   public:
      explicit AddAdditionalElecConvInfoProducer(const edm::ParameterSet&);
      ~AddAdditionalElecConvInfoProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      edm::InputTag eleSrc_;
      // ----------member data ---------------------------
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
AddAdditionalElecConvInfoProducer::AddAdditionalElecConvInfoProducer(const edm::ParameterSet& iConfig):
eleSrc_(iConfig.getParameter<edm::InputTag>("eleSrc"))
{
   //register your products

   produces<std::vector<pat::Electron> >();

}


AddAdditionalElecConvInfoProducer::~AddAdditionalElecConvInfoProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
AddAdditionalElecConvInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //B field
  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  float evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
  //electrons
  edm::Handle<std::vector<pat::Electron> > eles;
  iEvent.getByLabel(eleSrc_,eles);
  //Get the CTF tracks
  edm::Handle<reco::TrackCollection> tracks_h;
  iEvent.getByLabel("generalTracks", tracks_h);
  
  ConversionFinder convFinder;
  //returns the best candidate partner (see text below)
  ConversionInfo convInfo;
  std::vector<pat::Electron> newElecs;
  pat::Electron newElectron;
  for(std::vector<pat::Electron>::const_iterator it = eles->begin(); it != eles->end(); ++it){
   newElectron = *it;
   convInfo = convFinder.getConversionInfo(*it, tracks_h, evt_bField);
   newElectron.addUserFloat("deltaCotTheta",convInfo.dcot());
   newElectron.addUserFloat("deltaDistance",convInfo.dist());
   newElecs.push_back(newElectron);
 }

   // is put into the Event
   std::auto_ptr<std::vector<pat::Electron> > pOut(new std::vector<pat::Electron>(newElecs));
   iEvent.put(pOut);


/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
AddAdditionalElecConvInfoProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AddAdditionalElecConvInfoProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
AddAdditionalElecConvInfoProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
AddAdditionalElecConvInfoProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AddAdditionalElecConvInfoProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AddAdditionalElecConvInfoProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AddAdditionalElecConvInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddAdditionalElecConvInfoProducer);
