// -*- C++ -*-
//
// Package:    DoubleVecProducer
// Class:      DoubleVecProducer
// 
/**\class DoubleVecProducer DoubleVecProducer.cc MyProducers/DoubleVecProducer/src/DoubleVecProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Felix Hoehle
//         Created:  Thu Mar 29 00:39:36 CEST 2012
// $Id: DoubleVecProducer.cc,v 1.2 2012/08/11 16:34:50 fhohle Exp $
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

//
// class declaration
//

class DoubleVecProducer : public edm::EDProducer {
   public:
      explicit DoubleVecProducer(const edm::ParameterSet&);
      ~DoubleVecProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

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
DoubleVecProducer::DoubleVecProducer(const edm::ParameterSet& iConfig)
{
   //register your products
 produces<std::vector< double > >();

}


DoubleVecProducer::~DoubleVecProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DoubleVecProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   std::vector<double> weights;weights.push_back(1); 
   std::auto_ptr< std::vector< double > > pOut (new std::vector<double>(weights));
   edm::LogInfo("DoubleVecProducerPutEvt")<<"PU pOut weights "; for(unsigned int i = 0;i < pOut->size();i++){edm::LogInfo("DoubleVecProducerPutEvt")<<" "<<(*pOut)[i];} edm::LogInfo("DoubleVecProducerPutEvt")<<" end "<<std::endl;
   iEvent.put(pOut);


/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
DoubleVecProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DoubleVecProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
DoubleVecProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DoubleVecProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DoubleVecProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DoubleVecProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DoubleVecProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DoubleVecProducer);
