// -*- C++ -*-
//
// Package:    MyTTbarGenEvent10Parts
// Class:      MyTTbarGenEvent10Parts
// 
/**\class MyTTbarGenEvent10Parts MyTTbarGenEvent10Parts.cc MyProducers/MyTTbarGenEvent10Parts/src/MyTTbarGenEvent10Parts.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Felix Hoehle
//         Created:  Tue Mar 20 12:18:46 CET 2012
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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//
// class declaration
//

class MyTTbarGenEvent10Parts : public edm::EDProducer {
   public:
      explicit MyTTbarGenEvent10Parts(const edm::ParameterSet&);
      ~MyTTbarGenEvent10Parts();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

	const reco::GenParticle* firstParticleInChain(const reco::GenParticle* part);
	void ParticlePrint(const reco::GenParticle* part);
	std::vector<const reco::GenParticle*> ParticleChain(const reco::GenParticle* part);
	const reco::GenParticle* findDecayingParticle(std::vector<const reco::GenParticle*> parts);
	std::vector<reco::GenParticle> ConvertToObject(std::vector<const reco::GenParticle *> vecP);
      // ----------member data ---------------------------
};
std::vector<reco::GenParticle> MyTTbarGenEvent10Parts::ConvertToObject(std::vector<const reco::GenParticle *> vecP){
	std::vector<reco::GenParticle> vec;
	for(unsigned int i = 0; i < vecP.size(); ++i){
		vec.push_back(*vecP[i]);
	}
	return vec;
}
const reco::GenParticle* MyTTbarGenEvent10Parts::firstParticleInChain(const reco::GenParticle* part){
        if(part->mother() != 0) return (part->mother()->pdgId() == part->pdgId() ) ? NULL : part;
        else return part;
}
void MyTTbarGenEvent10Parts::ParticlePrint(const reco::GenParticle* part){

        edm::LogInfo  ("ParticleDecayProblem")<<"Printer Px: "<<part->px()<<"  Py: " <<part->py()<<"  Pz: "<<part->pz()<<"  E: "<<part->energy() <<"  pdgId "<<part->pdgId();
}
const reco::GenParticle* MyTTbarGenEvent10Parts::findDecayingParticle(std::vector<const reco::GenParticle*> parts){
	const reco::GenParticle* part = NULL;
	for(unsigned int i = 0; i < parts.size(); ++i){
	//	std::cout<<"num dau  "<<parts[i]->numberOfDaughters()<<std::endl;
                if(parts[i]->numberOfDaughters()>1){
                        if(part != NULL){ 
				edm::LogInfo  ("ParticleDecayProblem")<<"particle decays more than once!";
				ParticlePrint(parts[i]);
				ParticlePrint(part);
				edm::LogInfo  ("ParticleDecayProblem")<<"end print out"; }
                        part=parts[i];
                }
        }
	return part;
}
std::vector<const reco::GenParticle*> MyTTbarGenEvent10Parts::ParticleChain(const reco::GenParticle* part){
	std::vector<const reco::GenParticle*> particleChain;
	reco::GenParticleRefVector daughterRefs;
	particleChain.push_back(part);
	bool chainContinues=false;
	do{
		chainContinues=false;
		daughterRefs = particleChain[particleChain.size()-1]->daughterRefVector();
	        for(reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){
        	        if ( particleChain[particleChain.size()-1]->pdgId() == (*dau)->pdgId()){
	               	        chainContinues=true;
               	        	particleChain.push_back(&(**dau));
               	 	}
        	}
	}
	while(chainContinues);
	return particleChain;
}
//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MyTTbarGenEvent10Parts::MyTTbarGenEvent10Parts(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label*/
   produces<std::vector<std::vector<reco::GenParticle > > >("");
/* 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


MyTTbarGenEvent10Parts::~MyTTbarGenEvent10Parts()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MyTTbarGenEvent10Parts::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//  using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
*/
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel("genParticles",genParticles);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::vector<reco::GenParticle> temp;
   std::vector<const reco::GenParticle*> tops;
   std::vector<const reco::GenParticle*> antitops;
   std::auto_ptr<std::vector<std::vector< reco::GenParticle> > > pOut(new std::vector<std::vector< reco::GenParticle> >());
   unsigned int i = 0;
   for(reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++,i++){
	if(mcIter->pdgId() == 6 && firstParticleInChain(&(*mcIter)) > 0){
		//ParticlePrint(&(*mcIter));
		tops = MyTTbarGenEvent10Parts::ParticleChain(&(*mcIter));
   	}
   	if(mcIter->pdgId() == -6 && firstParticleInChain(&(*mcIter)) > 0){
		antitops = MyTTbarGenEvent10Parts::ParticleChain(&(*mcIter));
	}
   }
   //std::cout<<"tops "<<std::endl;
   //for(unsigned int i = 0;i<tops.size(); ++i){
   //	ParticlePrint(tops[i]);
   //}
   //std::cout<<"antitops "<<std::endl;
   //for(unsigned int i = 0;i<antitops.size(); ++i){
   //     ParticlePrint(antitops[i]);
   //}
   // tops okay?
   bool findThis=false;
   bool forAllOkay=true;
   if(tops.size() == 0 || antitops.size() == 0){ forAllOkay=false; std::cout<<"top size "<< tops.size() <<" antitop size "<< antitops.size() <<std::endl;}
   else{
   	for(unsigned int i = 0; i < tops[0]->numberOfMothers(); ++i){
		findThis=false;
		for(unsigned int j = 0; j < antitops[0]->numberOfMothers(); ++j){
			if( tops[0]->mother(i) == antitops[0]->mother(j)){ 
				findThis=true;
				break;
			}
		}
                if (!findThis) std::cout<<"mother check "<< i << "not found"<<std::endl;
		forAllOkay=forAllOkay&&findThis;
   	}
   }

   if(forAllOkay){
	// find decaying tops and their decay products
	const reco::GenParticle* decayingTop = findDecayingParticle(tops);
	//ParticlePrint(decayingTop);
	const reco::GenParticle* b = NULL;
	const reco::GenParticle* Wplus = NULL; 
	reco::GenParticleRefVector topdaughterRefs = decayingTop->daughterRefVector();
	for(reco::GenParticleRefVector::const_iterator dau = topdaughterRefs.begin(); dau!=topdaughterRefs.end(); ++dau){
		if((*dau)->pdgId() == 24){
			Wplus =dynamic_cast <const reco::GenParticle *>(&(**dau));
		}
                if((*dau)->pdgId() <= 6 && (*dau)->pdgId() >= 1){
                        b = dynamic_cast <const reco::GenParticle *>(&(**dau));
                }
	}
	//ParticlePrint(b);
	//ParticlePrint(Wplus);
	const reco::GenParticle* decayingAntitop = findDecayingParticle(antitops);
        const reco::GenParticle* bBar = NULL;
        const reco::GenParticle* Wminus = NULL;
        //ParticlePrint(decayingAntitop);
	reco::GenParticleRefVector antitopdaughterRefs = decayingAntitop->daughterRefVector();
        for(reco::GenParticleRefVector::const_iterator dau = antitopdaughterRefs.begin(); dau!=antitopdaughterRefs.end(); ++dau){
                if((*dau)->pdgId() == -24){
                        Wminus = dynamic_cast <const reco::GenParticle *>(&(**dau));
                }
                if((*dau)->pdgId() >= -6 && (*dau)->pdgId() <= -1){
                        bBar = dynamic_cast <const reco::GenParticle *>(&(**dau));
                }
        }
	//ParticlePrint(Wminus);
	//ParticlePrint(bBar);
	// topDecaysOkay ?
	if(b!=0 && bBar !=0 && Wplus != 0 && Wminus != 0){
		std::vector<const reco::GenParticle*> bBars =ParticleChain(bBar); 
	   	std::vector<const reco::GenParticle*> Wminuss = ParticleChain(Wminus);
	   	std::vector<const reco::GenParticle*> bs = ParticleChain(b);
	   	std::vector<const reco::GenParticle*> Wpluss = ParticleChain(Wplus);
		//std::cout<<"bBars "<<std::endl;
		//for (unsigned int i = 0; i < bBars.size(); ++i){
		//	std::cout<<"bB ";
		//	ParticlePrint(bBars[i]);
		//}	
	
		// find decaying Ws
		const reco::GenParticle * decayingWplus = findDecayingParticle(Wpluss);
		const reco::GenParticle * WplusD1 = NULL;
		const reco::GenParticle * WplusD2 = NULL;
		reco::GenParticleRefVector WplusdaughterRefs = decayingWplus->daughterRefVector();
	        for(reco::GenParticleRefVector::const_iterator dau = WplusdaughterRefs.begin(); dau !=WplusdaughterRefs.end();++dau ){
			int pdgId = (*dau)->pdgId();
	                if(abs(pdgId) <= 6 || (abs(pdgId) >= 11 && abs(pdgId) <= 16)){
				if(pdgId < 0) WplusD1 = dynamic_cast <const reco::GenParticle *>(&(**dau));
				if(pdgId > 0) WplusD2 = dynamic_cast <const reco::GenParticle *>(&(**dau));
	                }
	        }
		//ParticlePrint(WplusD1);
		std::vector<const reco::GenParticle*> WplusD1s = ParticleChain(WplusD1);
		//ParticlePrint(WplusD2);
		std::vector<const reco::GenParticle*> WplusD2s = ParticleChain(WplusD2);
                //for (unsigned int i = 0; i < WplusD1s.size(); ++i){
                //        std::cout<<"WplusD1s ";
                //        ParticlePrint(WplusD1s[i]);
                //} 
                //for (unsigned int i = 0; i < WplusD2s.size(); ++i){
                //        std::cout<<"WplusD2s ";
                //        ParticlePrint(WplusD2s[i]);
                //}

	        const reco::GenParticle * decayingWminus = findDecayingParticle(Wminuss);
	        const reco::GenParticle * WminusD1 = NULL;
	        const reco::GenParticle * WminusD2 = NULL;
		reco::GenParticleRefVector WminusdaughterRefs = decayingWminus->daughterRefVector();
                for(reco::GenParticleRefVector::const_iterator dau = WminusdaughterRefs.begin(); dau != WminusdaughterRefs.end(); ++dau){
	                int pdgId = (*dau)->pdgId();
	                if(abs(pdgId) <= 6 || (abs(pdgId) >= 11 && abs(pdgId) <= 16)){
	                        if(pdgId < 0) WminusD1 = dynamic_cast <const reco::GenParticle *>(&(**dau));
                        	if(pdgId > 0) WminusD2 = dynamic_cast <const reco::GenParticle *>(&(**dau));
	            	}
		}
		//ParticlePrint(WminusD1);
		std::vector<const reco::GenParticle*> WminusD1s = ParticleChain(WminusD1);
		//ParticlePrint(WminusD2);
		std::vector<const reco::GenParticle*> WminusD2s = ParticleChain(WminusD2);
   		pOut->push_back(ConvertToObject(tops));
   		pOut->push_back(ConvertToObject(bs));
   		pOut->push_back(ConvertToObject(Wpluss));
   		pOut->push_back(ConvertToObject(WplusD1s));
   		pOut->push_back(ConvertToObject(WplusD2s));
   		//
   		pOut->push_back(ConvertToObject(antitops));
   		pOut->push_back(ConvertToObject(bBars));
   		pOut->push_back(ConvertToObject(Wminuss));
   		pOut->push_back(ConvertToObject(WminusD1s));
   		pOut->push_back(ConvertToObject(WminusD2s));

        } 
   }
   else std::cout<<"not found: tops are not okay"<<std::endl;
   iEvent.put(pOut);

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
MyTTbarGenEvent10Parts::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyTTbarGenEvent10Parts::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
MyTTbarGenEvent10Parts::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MyTTbarGenEvent10Parts::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyTTbarGenEvent10Parts::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyTTbarGenEvent10Parts::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyTTbarGenEvent10Parts::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyTTbarGenEvent10Parts);
