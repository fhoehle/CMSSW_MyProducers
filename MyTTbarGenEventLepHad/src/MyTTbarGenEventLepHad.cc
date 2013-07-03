// -*- C++ -*-
//
// Package:    MyTTbarGenEventLepHad
// Class:      MyTTbarGenEventLepHad
// 
/**\class MyTTbarGenEventLepHad MyTTbarGenEventLepHad.cc MyProducers/MyTTbarGenEventLepHad/src/MyTTbarGenEventLepHad.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Felix Hoehle
//         Created:  Wed Jan 11 11:29:35 CET 2012
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

class MyTTbarGenEventLepHad : public edm::EDProducer {
   public:
      explicit MyTTbarGenEventLepHad(const edm::ParameterSet&);
      ~MyTTbarGenEventLepHad();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag genTag_;

	const reco::GenParticle* firstParticleInChain(const reco::GenParticle* part);
	std::vector<std::vector<const reco::GenParticle*> > isrParticles(const reco::GenParticle* part);
	void ParticlePrint(const reco::GenParticle* part);
	std::vector<std::vector< reco::GenParticle > > PointerToCollection(std::vector<std::vector<const reco::GenParticle*> > PointerCollection);
	bool DecayingW (reco::GenParticle* part);
	void FindRightDecayProds (reco::GenParticle* part, reco::GenParticle &part1, reco::GenParticle &part2);
	edm::InputTag genPartTag_;

      // ----------member data ---------------------------
};
bool MyTTbarGenEventLepHad::DecayingW (reco::GenParticle* part){
	reco::GenParticleRefVector daughterRefs;
	daughterRefs = part->daughterRefVector();
	int numDau = part->numberOfDaughters();
	if(numDau <= 0) return false;
	else if (abs(part->pdgId()) != 24 ) return false;
	else{
		//std::cout<<"DecayingW ";
		//ParticlePrint(part);
		//std::cout<<"DecayingW daughters"<<part->numberOfDaughters() <<std::endl;

		int counterQuarkDaughters = 0;
		for(reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){				
			//ParticlePrint(&(**dau));
			//std::cout<<"DecayingW daughters pdgId "<<abs((*dau)->pdgId())<<std::endl;
			if( (abs((*dau)->pdgId()) >= 1 && abs((*dau)->pdgId()) <= 6) ||  (abs((*dau)->pdgId()) >= 11 && abs((*dau)->pdgId()) <= 16)  ){ counterQuarkDaughters++;}
		}
		//std::cout<<"counterQuarkDaughters "<<counterQuarkDaughters <<std::endl;
		//std::cout<<"DecayingW finish"<<part->numberOfDaughters() <<std::endl;
		if(counterQuarkDaughters == 2) return true;
		else return false;
	}

}
void MyTTbarGenEventLepHad::FindRightDecayProds (reco::GenParticle* part, reco::GenParticle &part1, reco::GenParticle &part2){
	if( DecayingW(part) ){
		reco::GenParticleRefVector daughterRefs;
		daughterRefs = part->daughterRefVector();
		int flag = 0;
		for(reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){
			if(abs((*dau)->pdgId()) != 24 && ( (abs((*dau)->pdgId()) <= 6 && abs((*dau)->pdgId()) >= 1 ) || (abs((*dau)->pdgId()) >= 11 && abs((*dau)->pdgId()) <= 16) ) && flag == 1 ){
				part1	= **(dau);
				flag++;
			}
			if(abs((*dau)->pdgId()) != 24 && ( (abs((*dau)->pdgId()) <= 6 && abs((*dau)->pdgId()) >= 1) || (abs((*dau)->pdgId()) >= 11 && abs((*dau)->pdgId()) <= 16) ) && flag == 0 ){
				part2	= **(dau);
				flag++;
			} 
		}
		if (flag != 2) std::cout<<"Very strange W Decay"<<std::endl;
	}
	else std::cout<<"Wrong usage FindRightDecayProds"<<std::endl;

}
const reco::GenParticle* MyTTbarGenEventLepHad::firstParticleInChain(const reco::GenParticle* part){
	if(part->mother() != 0) return (part->mother()->pdgId() == part->pdgId() ) ? NULL : part;
	else return part;
}
void MyTTbarGenEventLepHad::ParticlePrint(const reco::GenParticle* part){

	std::cout<<"Printer Px: "<<part->px()<<"  Py: " <<part->py()<<"  Pz: "<<part->pz()<<"  E: "<<part->energy() <<"  pdgId "<<part->pdgId() <<std::endl;
	}
std::vector<std::vector< reco::GenParticle > > MyTTbarGenEventLepHad::PointerToCollection(std::vector<std::vector<const reco::GenParticle*> > PointerCollection){
	std::vector<std::vector< reco::GenParticle > > collection;
	std::vector<reco::GenParticle > tempCollection;
	for(uint i =0; i <PointerCollection.size(); i++ ){
		tempCollection.clear();
		for(uint j =0; j <PointerCollection[i].size(); j++ ){
			tempCollection.push_back(*(PointerCollection[i][j]));
		}
		collection.push_back(tempCollection);
	}
	return collection;
 }
	
std::vector<std::vector<const reco::GenParticle*> > MyTTbarGenEventLepHad::isrParticles(const reco::GenParticle* part){
	std::vector<std::vector<const reco::GenParticle*> > result;
	std::vector<const reco::GenParticle*> temp;
	reco::GenParticleRefVector daughterRefs;
	const reco::GenParticle* temppart = part; 
	const reco::GenParticle* nexttempPart;

	for(unsigned int i = 0;i < part->numberOfMothers(); i++ ){
		if(part->pdgId() == (part->mother(i))->pdgId()){std::cout<<"Wrong usage isrParticles() "<<std::endl; return result;}
	}

	temp.push_back(temppart);
	result.push_back(temp);

	int idx = -1;
	int z = 0;
	do{	
		temp.resize(0);
		idx = -1; z = 0; nexttempPart = NULL;
		daughterRefs = temppart->daughterRefVector();
		for( reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau,z++){
			if(temppart->pdgId() == (*dau)->pdgId()){idx = z; nexttempPart = &(**dau);} 
		}
		if(idx != -1){
			temp.push_back(nexttempPart);
			z = 0;
			for( reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau,z++){
				if(z != idx){temp.push_back((&**dau)); }
			}
			temppart = nexttempPart;
			result.push_back(temp);
		}
	}
	while( idx != -1);

	return result;
}

//
// constants, enums and typedefs



//
// static data member definitions
//
//
// constructors and destructor
//
MyTTbarGenEventLepHad::MyTTbarGenEventLepHad(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
	genTag_ = iConfig.getUntrackedParameter<edm::InputTag>("genTag");  
	produces<std::vector<reco::GenParticle> >().setBranchAlias( "MyTTbarGenEvent");
}


MyTTbarGenEventLepHad::~MyTTbarGenEventLepHad()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MyTTbarGenEventLepHad::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	edm::Handle< edm::View<reco::GenParticle> > genParts;
    	iEvent.getByLabel(genTag_, genParts);



	std::auto_ptr< std::vector< std::vector< std::vector<reco::GenParticle> > > >  MyTTbarEvent ( new std::vector<std::vector< std::vector<reco::GenParticle> > > );

	int NGenParts = genParts->size();
	reco::GenParticle tempgenPart;
	std::vector<std::vector<const reco::GenParticle*> > topChainPointer;
	std::vector<std::vector<reco::GenParticle> > topChain;
	std::vector<reco::GenParticle> topCollection;
	reco::GenParticle top;

	std::vector<std::vector<const reco::GenParticle*> > antitopChainPointer;
	std::vector<std::vector<reco::GenParticle> > antitopChain;
	std::vector<reco::GenParticle> antitopCollection;
	reco::GenParticle antitop;

	std::vector<std::vector<reco::GenParticle> > WminusChain;
	reco::GenParticle Wminus;
	std::vector<std::vector<reco::GenParticle> > WplusChain;
	reco::GenParticle Wplus;
	std::vector<std::vector<reco::GenParticle> > bChain;
	reco::GenParticle b;
	std::vector<std::vector<reco::GenParticle> > bbarChain;
	reco::GenParticle bbar;

	reco::GenParticle WplusDecay1;
	std::vector<std::vector<reco::GenParticle> > WplusDecay1Chain;
	reco::GenParticle WplusDecay2;
	std::vector<std::vector<reco::GenParticle> > WplusDecay2Chain;

	reco::GenParticle WminusDecay1;
	std::vector<std::vector<reco::GenParticle> > WminusDecay1Chain;
	reco::GenParticle WminusDecay2;
	std::vector<std::vector<reco::GenParticle> > WminusDecay2Chain;


	std::vector<std::vector<const reco::GenParticle*> > tempChainPointer;
	std::vector<reco::GenParticle> tempCollection;

	int NumTops = 0;
	for(int idx=0; idx<NGenParts; ++idx) {
		tempgenPart =  (*genParts)[idx];
		if(abs(tempgenPart.pdgId()) == 6 && tempgenPart.mother() != 0 &&  (tempgenPart.mother())->pdgId() != tempgenPart.pdgId() ){
			if(tempgenPart.pdgId() == 6) top = tempgenPart;
			else if(tempgenPart.pdgId() == -6) antitop = tempgenPart;
			else std::cout<<"This should not happen!"<<std::endl;
			NumTops++;
		}
		
	}
	if(NumTops >= 2){
		topChainPointer = isrParticles(&(top));
		topChain = PointerToCollection(topChainPointer);

		antitopChainPointer = isrParticles(&(antitop));
		antitopChain = PointerToCollection(antitopChainPointer);
		// Wminus search
		reco::GenParticle lastPart;
		for(uint i =0; i < antitopChain[antitopChain.size()-1].size(); i++ ){
			if(antitopChain[antitopChain.size()-1][i].pdgId() == antitopChain[0][0].pdgId()){
				lastPart = antitopChain[antitopChain.size()-1][i];
			}
		}
		//std::cout<<"lastPart: antitopChain"<<std::endl;
		//ParticlePrint(&lastPart);
		if(lastPart.pdgId() == antitopChain[0][0].pdgId()){
			reco::GenParticleRefVector daughterRefs;
			daughterRefs = lastPart.daughterRefVector();
			for( reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){
				if((*dau)->pdgId() == -24 ){
					Wminus	= **(dau);
					//ParticlePrint(&Wminus);
				}
			}
		}
		// Wplus search
		reco::GenParticle lastPart1;
		for(uint i =0; i < topChain[topChain.size()-1].size(); i++ ){
			if(topChain[topChain.size()-1][i].pdgId() == topChain[0][0].pdgId()){
				lastPart1 = topChain[topChain.size()-1][i];
			}
		}
		//std::cout<<"lastPart1: topChain"<<std::endl;
		//ParticlePrint(&lastPart1);
		if(lastPart1.pdgId() == topChain[0][0].pdgId()){
			reco::GenParticleRefVector daughterRefs;
			daughterRefs = lastPart1.daughterRefVector();
			for( reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){
				if((*dau)->pdgId() == 24 ){
					Wplus	= **(dau);
					//ParticlePrint(&Wplus);
				}
			}
		}
		// b search
		reco::GenParticle lastPart1b;
		for(uint i =0; i < topChain[topChain.size()-1].size(); i++ ){
			if(topChain[topChain.size()-1][i].pdgId() == topChain[0][0].pdgId()){
				lastPart1b = topChain[topChain.size()-1][i];
				
			}
		}
		//std::cout<<"lastPart1b: topChain"<<std::endl;
		//ParticlePrint(&lastPart1b);
		//std::cout<<"lastPart1b: topChaindatwillich"<<std::endl;
		if(lastPart1b.pdgId() == topChain[0][0].pdgId()){
			reco::GenParticleRefVector daughterRefs;
			daughterRefs = lastPart1b.daughterRefVector();
			for( reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){
				if((*dau)->pdgId() == 5 || (*dau)->pdgId() == 3 || (*dau)->pdgId() == 1){
					b	= **(dau);
					//ParticlePrint(&b);
				}
			}
		}
		// bbar search
		reco::GenParticle lastPartbbar;
		for(uint i =0; i < antitopChain[antitopChain.size()-1].size(); i++ ){
			if(antitopChain[antitopChain.size()-1][i].pdgId() == antitopChain[0][0].pdgId()){
				lastPartbbar = antitopChain[antitopChain.size()-1][i];
			}
		}
		//std::cout<<"lastPartbbar: antitopChain"<<std::endl;
		//ParticlePrint(&lastPartbbar);
		if(lastPartbbar.pdgId() == antitopChain[0][0].pdgId()){
			reco::GenParticleRefVector daughterRefs;
			daughterRefs = lastPartbbar.daughterRefVector();
			for( reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){
				if((*dau)->pdgId() == -5  || (*dau)->pdgId() == -3 || (*dau)->pdgId() == -1){
					bbar	= **(dau);
					//ParticlePrint(&bbar);
				}
			}
		}

	}
	else std::cout<<"Not enough Tops found"<<std::endl;
	//std::cout<<"troet"<<std::endl;


	if(Wminus.pdgId() != -24) std::cout<<"huch Wminus "<< std::endl;
	else{	
		tempChainPointer = isrParticles(&(Wminus));
		WminusChain = PointerToCollection(tempChainPointer);
		int numberOfDecayingW = 0;
		for(uint i =0; i < WminusChain.size(); i++){
			if( DecayingW( &(WminusChain[i][0])) ){
			FindRightDecayProds( &(WminusChain[i][0]) ,WminusDecay1,WminusDecay2);
			numberOfDecayingW++;
			}
		}
		if(numberOfDecayingW >=2 ) std::cout<<"Too many decaying W minus"<<std::endl;

// 		reco::GenParticleRefVector daughterRefs;
// 		daughterRefs = WminusChain[WminusChain.size()-1][0].daughterRefVector();
// 		int flag = 0;
		//std::cout<<"WminusDecay"<<std::endl;
// 		for( reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){
// 			if((*dau)->pdgId() != -24 && flag == 1 ){
// 				WminusDecay2	= **(dau);
// 				//ParticlePrint(&WminusDecay2);
// 				flag++;
// 			}
// 			if((*dau)->pdgId() != -24 && flag == 0 ){
// 				WminusDecay1	= **(dau);
// 				//ParticlePrint(&WminusDecay1);
// 				flag++;
// 			}
// 
// 		}
		//ParticlePrint(&WminusDecay1);
		//ParticlePrint(&WminusDecay2);
		
	}
	//std::cout<<"troet11"<<std::endl;
	if(WminusDecay1.pdgId() != 0) WminusDecay1Chain = PointerToCollection(isrParticles(&(WminusDecay1)));
	//else std::cout<<"Holla wasn los WminusDecay1"<<std::endl;

	//std::cout<<"troet11a"<<std::endl;
	if(WminusDecay2.pdgId() != 0) WminusDecay2Chain = PointerToCollection(isrParticles(&(WminusDecay2)));
	//else std::cout<<"Holla wasn los WminusDecay2"<<std::endl;

	//std::cout<<"troet111"<<std::endl;
	if(Wplus.pdgId() != 24) std::cout<<"huch Wplus "<< std::endl;
	else{	
		tempChainPointer = isrParticles(&(Wplus));
		WplusChain = PointerToCollection(tempChainPointer);

		int numberOfDecayingW = 0;
		for(uint i =0; i < WplusChain.size(); i++){
			if(DecayingW( &(WplusChain[i][0])) ){
				FindRightDecayProds( &(WplusChain[i][0]),WplusDecay1,WplusDecay2);
				numberOfDecayingW++;
				//std::cout<<"Decaying W plus"<<std::endl;
				//ParticlePrint(&(WplusChain[i][0]));
			}
			else{
				//std::cout<<"This W plus not "<<WplusChain[i][0].numberOfDaughters() <<std::endl;
				//ParticlePrint(&(WplusChain[i][0]));	

			}
		}
		if(numberOfDecayingW >=2 || numberOfDecayingW == 0){ 
			std::cout<<"Too many decaying W plus or no/not enough W plus Decay: "<<numberOfDecayingW<<std::endl;
			ParticlePrint(&Wplus);
		}

// 		reco::GenParticleRefVector daughterRefs;
// 		daughterRefs = WplusChain[WplusChain.size()-1][0].daughterRefVector();
// 		int flag = 0;
// 		for( reco::GenParticleRefVector::const_iterator dau = daughterRefs.begin(); dau!=daughterRefs.end(); ++dau){
// 			if((*dau)->pdgId() != 24 && flag == 1 ){
// 				WplusDecay2	= **(dau);
// 				flag++;
// 			}
// 			if((*dau)->pdgId() != 24 && flag == 0 ){
// 				WplusDecay1	= **(dau);
// 				flag++;
// 			}
// 
// 		}
		//ParticlePrint(&WplusDecay1);
		//ParticlePrint(&WplusDecay2);
	}
	if(WplusDecay1.pdgId() != 0) WplusDecay1Chain = PointerToCollection(isrParticles(&(WplusDecay1)));
	//else std::cout<<"Holla wasn los WplusDecay1"<<std::endl;
	if(WplusDecay2.pdgId() != 0) WplusDecay2Chain = PointerToCollection(isrParticles(&(WplusDecay2)));
	//else std::cout<<"Holla wasn los WplusDecay2"<<std::endl;
	//std::cout<<"troet1ee"<<std::endl;
	if(bbar.pdgId() != -5 && bbar.pdgId() != -3 && bbar.pdgId() != -1) std::cout<<"huch bbar "<< std::endl;
	else{	
		tempChainPointer = isrParticles(&(bbar));
		bbarChain = PointerToCollection(tempChainPointer);
	}
	//std::cout<<"huch b wasn los "<< std::endl;
	//ParticlePrint(&b);
	if(b.pdgId() != 5 && b.pdgId() != 3 && b.pdgId() != 1) std::cout<<"huch b "<< std::endl;
	else{	
		tempChainPointer = isrParticles(&(b));
		bChain = PointerToCollection(tempChainPointer);
	}

	//std::cout<<"Num Tops: "<<NumTops<<std::endl;
	(*MyTTbarEvent).push_back(topChain);
	(*MyTTbarEvent).push_back(bChain);
	(*MyTTbarEvent).push_back(WplusChain);
	(*MyTTbarEvent).push_back(WplusDecay1Chain);
	(*MyTTbarEvent).push_back(WplusDecay2Chain);
	//
	(*MyTTbarEvent).push_back(antitopChain);
	(*MyTTbarEvent).push_back(bbarChain);
	(*MyTTbarEvent).push_back(WminusChain);
	(*MyTTbarEvent).push_back(WminusDecay1Chain);
	(*MyTTbarEvent).push_back(WminusDecay2Chain);
	//ParticlePrint(&topChain[0][0]);
	//ParticlePrint(&bChain[0][0]);
	//ParticlePrint(&WplusChain[0][0]);
	//ParticlePrint(&WplusDecay1Chain[0][0]);
	//ParticlePrint(&WplusDecay2Chain[0][0]);
	//ParticlePrint(&antitopChain[0][0]);
	//ParticlePrint(&bbarChain[0][0]);
	//ParticlePrint(&WminusChain[0][0]);
	//ParticlePrint(&WminusDecay1Chain[0][0]);
	//ParticlePrint(&WminusDecay2Chain[0][0]);	
	
	bool topleptonic = (fabs(WplusDecay1Chain[0][0].pdgId()) > 10);
	//std::cout<<"top pdgId() "<<topChain[0][0].pdgId()<<" W "<< WplusChain[0][0].pdgId()<<" b "<<bChain[0][0].pdgId()<<" Wdecays "<<WplusDecay1Chain[0][0].pdgId()<<" "<<WplusDecay2Chain[0][0].pdgId()<<std::endl;
	//std::cout<<"antitop pdgId() "<<antitopChain[0][0].pdgId()<<" W "<< WminusChain[0][0].pdgId()<<" b "<<bbarChain[0][0].pdgId()<<" Wdecays "<<WminusDecay1Chain[0][0].pdgId()<<" "<<WminusDecay2Chain[0][0].pdgId()<<std::endl;

	reco::GenParticle lepton = ( (topleptonic) ? (WplusDecay1Chain[0][0].pdgId()%2 == 1 ? WplusDecay1Chain[0][0]: WplusDecay2Chain[0][0]) :(WminusDecay1Chain[0][0].pdgId()%2 == 1 ? WminusDecay1Chain[0][0] : WminusDecay2Chain[0][0])  );
	reco::GenParticle neutrino = ( (topleptonic) ? (WplusDecay1Chain[0][0].pdgId()%2 == 1 ? WplusDecay2Chain[0][0]: WplusDecay1Chain[0][0]) :(WminusDecay1Chain[0][0].pdgId()%2 == 1 ? WminusDecay2Chain[0][0] : WminusDecay1Chain[0][0])  );	
	reco::GenParticle ulikequark = ( (topleptonic) ?  (WminusDecay1Chain[0][0].pdgId()%2 == 1 ? WminusDecay2Chain[0][0] : WminusDecay1Chain[0][0]) : (WplusDecay1Chain[0][0].pdgId()%2 == 1 ? WplusDecay2Chain[0][0]: WplusDecay1Chain[0][0])  );
	reco::GenParticle dlikequark = ( (topleptonic) ?  (WminusDecay1Chain[0][0].pdgId()%2 == 1 ? WminusDecay1Chain[0][0] : WminusDecay2Chain[0][0]) : (WplusDecay1Chain[0][0].pdgId()%2 == 1 ? WplusDecay1Chain[0][0]: WplusDecay2Chain[0][0])  );

	reco::GenParticle blep = topleptonic ? bChain[0][0] : bbarChain[0][0];
	reco::GenParticle bhad = topleptonic ? bbarChain[0][0] : bChain[0][0]; 
	reco::GenParticle Wlep = topleptonic ? WplusChain[0][0] : WminusChain[0][0];
	reco::GenParticle Whad = topleptonic ?  WminusChain[0][0] : WplusChain[0][0];
	reco::GenParticle toplep = topleptonic ? topChain[topChain.size()-1][0] : antitopChain[antitopChain.size()-1][0];
	reco::GenParticle tophad = topleptonic ? antitopChain[antitopChain.size()-1][0] : topChain[topChain.size()-1][0];

std::auto_ptr< std::vector<reco::GenParticle> >  MyTTbarEventVec ( new std::vector<reco::GenParticle> );
	// std::cout<<"topleptonic "<<topleptonic<<std::endl;
        //std::cout<<"top toplep "<<toplep.pdgId()<<" W "<< Wlep.pdgId()<<" b "<<blep.pdgId()<<" Wdecays "<<lepton.pdgId()<<" "<<neutrino.pdgId()<<std::endl;
        //std::cout<<"top had "<<tophad.pdgId()<<" W "<< Whad.pdgId()<<" b "<<bhad.pdgId()<<" Wdecays "<<dlikequark.pdgId()<<" "<<ulikequark.pdgId()<<std::endl;


	(*MyTTbarEventVec).push_back(toplep);
	(*MyTTbarEventVec).push_back(blep);
	(*MyTTbarEventVec).push_back(Wlep);
	(*MyTTbarEventVec).push_back(lepton);
	(*MyTTbarEventVec).push_back(neutrino);

	(*MyTTbarEventVec).push_back(tophad),
	(*MyTTbarEventVec).push_back(bhad);
	(*MyTTbarEventVec).push_back(Whad);
	(*MyTTbarEventVec).push_back(dlikequark);
	(*MyTTbarEventVec).push_back(ulikequark);
//	(*MyTTbarEvent).push_back(&topChain);
//	(*MyTTbarEvent).push_back(bChain);
//	(*MyTTbarEvent).push_back(WplusChain);
//	(*MyTTbarEvent).push_back(WplusDecay1Chain);
//	(*MyTTbarEvent).push_back(WplusDecay2Chain);
	//
//	(*MyTTbarEvent).push_back(antitopChain);
//	(*MyTTbarEvent).push_back(bbarChain);
//	(*MyTTbarEvent).push_back(WminusChain);
//	(*MyTTbarEvent).push_back(WminusDecay1Chain);
//	(*MyTTbarEvent).push_back(WminusDecay2Chain);
//	(*MyTTbarEventVec)


	iEvent.put(MyTTbarEventVec);// MyTTbarEvent antitopChain
}

// ------------ method called once each job just before starting event loop  ------------
void 
MyTTbarGenEventLepHad::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyTTbarGenEventLepHad::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
MyTTbarGenEventLepHad::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MyTTbarGenEventLepHad::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyTTbarGenEventLepHad::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyTTbarGenEventLepHad::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyTTbarGenEventLepHad::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyTTbarGenEventLepHad);
