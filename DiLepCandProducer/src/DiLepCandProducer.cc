// -*- C++ -*-
//
// Package:    DiLepCandProducer
// Class:      DiLepCandProducer
// 
/**\class DiLepCandProducer DiLepCandProducer.cc CMSSW_MyProducers/DiLepCandProducer/src/DiLepCandProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Felix Hoehle
//         Created:  Thu Oct  4 23:26:24 CEST 2012
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

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "PhysicsTools/PatUtils/interface/PATDiObjectProxy.h"
#include "DataFormats/Math/interface/deltaR.h"
//
// class declaration
//

class DiLepCandProducer : public edm::EDProducer {
   public:
      explicit DiLepCandProducer(const edm::ParameterSet&);
      ~DiLepCandProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag srcColl1_, srcColl2_;
	std::string cut_;
	StringCutObjectSelector< pat::DiObjectProxy > 	pairCut_;
//	bool HighestPtSum( pat::CompositeCandidate p1, pat::CompositeCandidate p2 );
      // ----------member data ---------------------------
};
struct HighestPtSum {
 bool operator()( pat::CompositeCandidate p1, pat::CompositeCandidate p2 ) const {
//        std::cout<<"ptsums "<< (&p1)->daughter("l1")->p4().Pt()+(&p1)->daughter("l2")->p4().Pt() << "  " <<  (&p2)->daughter("l1")->p4().Pt()+(&p2)->daughter("l2")->p4().Pt() <<std::endl;
        return (&p1)->daughter("p1")->p4().Pt()+(&p1)->daughter("p2")->p4().Pt() >  (&p2)->daughter("p1")->p4().Pt()+(&p2)->daughter("p2")->p4().Pt() ;
    }
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
DiLepCandProducer::DiLepCandProducer(const edm::ParameterSet& iConfig):
srcColl1_(iConfig.getParameter<edm::InputTag>("srcColl1")),
srcColl2_(iConfig.getParameter<edm::InputTag>("srcColl2")),
cut_(iConfig.getParameter<std::string>("cut")),
pairCut_(iConfig.getParameter<std::string>("pairCut"))
{
   produces<std::vector<pat::CompositeCandidate> >();
}


DiLepCandProducer::~DiLepCandProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DiLepCandProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
  edm::Handle<edm::View<reco::Candidate> > coll1;
  iEvent.getByLabel(srcColl1_,coll1);
  edm::Handle<edm::View<reco::Candidate> > coll2;
  iEvent.getByLabel(srcColl2_,coll2);
  std::vector<pat::CompositeCandidate> diLepCands;
//  std::cout<<"colls "<<coll1->size() <<"  "<<coll2->size()<<std::endl;
  StringCutObjectSelector<pat::CompositeCandidate> select( cut_.c_str(),true);
  for(edm::View<reco::Candidate>::const_iterator itPart1= coll1->begin(); itPart1 != coll1->end(); ++itPart1){
    for(edm::View<reco::Candidate>::const_iterator itPart2= coll2->begin(); itPart2 != coll2->end(); ++itPart2){
      if(pairCut_(pat::DiObjectProxy(*itPart1,*itPart2))){
	pat::CompositeCandidate diLepCand;
	diLepCand.addDaughter( *itPart1, "p1");
    	diLepCand.addDaughter( *itPart2, "p2");
   	diLepCand.setP4(itPart1->p4()+itPart2->p4());
   	AddFourMomenta addp4;
//   addp4.set( diLepCand );
	if(select(diLepCand)) diLepCands.push_back(diLepCand);
	 }
	}
  }
//   std::cout<<" size cand Coll "<<diLepCands.size()<<std::endl;
  std::sort( diLepCands.begin(),  diLepCands.end(),HighestPtSum());
//  std::cout<<"donesort"<<std::endl;
  std::auto_ptr<std::vector<pat::CompositeCandidate> > pOut( new std::vector<pat::CompositeCandidate>());
  if(diLepCands.size()> 0) pOut->push_back(diLepCands[0]);
  iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void 
DiLepCandProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLepCandProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
DiLepCandProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DiLepCandProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DiLepCandProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DiLepCandProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiLepCandProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiLepCandProducer);
