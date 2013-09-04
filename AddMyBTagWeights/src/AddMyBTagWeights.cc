// -*- C++ -*-
//
// Package:    AddMyBTagWeights
// Class:      AddMyBTagWeights
// 
/**\class AddMyBTagWeights AddMyBTagWeights.cc CMSSW_MyProducers/AddMyBTagWeights/src/AddMyBTagWeights.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Felix Hoehle
//         Created:  Fri Nov  9 21:55:42 CET 2012
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
//
// class declaration
//

class AddMyBTagWeights : public edm::EDProducer {
   public:
      explicit AddMyBTagWeights(const edm::ParameterSet&);
      ~AddMyBTagWeights();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  edm::InputTag jetSrc_;
  TFile * histFile;
  TFile * btagEffBTagFile;
  TH2D * btagHist;
  TH2F * h_eff_cq; TH2F * h_eff_lq ; TH2F * h_eff_bq ;

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
AddMyBTagWeights::AddMyBTagWeights(const edm::ParameterSet& iConfig):
jetSrc_(iConfig.getParameter<edm::InputTag>("jetSrc"))
{
   //register your products
	histFile = new TFile((std::string(getenv("CMSSW_BASE"))+std::string("/src/CMSSW_MyProducers/AddMyBTagWeights/data/efficacite_btag.root")).c_str());btagHist =  new TH2D("btagHist","",3,1,2,3,1,2);
  histFile->GetObject("h_MISTAGCSVL_BTAGLEFFCORR",btagHist);
	btagEffBTagFile = new TFile((std::string(getenv("CMSSW_BASE"))+std::string("/src/CMSSW_MyProducers/AddMyBTagWeights/data/eff_from_ttmadgraph_summer11_multiwp.root")).c_str());
  h_eff_cq = new TH2F("h_eff_cq","",1,2,3,1,2,3); h_eff_lq = new TH2F("h_eff_lq","",1,2,3,1,2,3); h_eff_bq  = new TH2F("h_eff_bq","",1,2,3,1,2,3);
  btagEffBTagFile->GetObject("algo_6_discri_0.244/h_eff_cq;1",h_eff_cq); btagEffBTagFile->GetObject("algo_6_discri_0.244/h_eff_bq;1",h_eff_bq); btagEffBTagFile->GetObject("algo_6_discri_0.244/h_eff_lq;1",h_eff_lq);
   produces<std::vector<pat::Jet> >();
/*
   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


AddMyBTagWeights::~AddMyBTagWeights()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
AddMyBTagWeights::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   // required number of BJets
   int n_bjets_  = 2;

   std::vector<double> btagWeights;std::vector<double> proba_jetb;
   std::vector<double> btagEff;

   edm::Handle<std::vector<pat::Jet> > myJets;
   iEvent.getByLabel(jetSrc_,myJets);
   std::vector<pat::Jet> newJets; pat::Jet newJet;
   int jetFlav; double jetPt; double jetEta; double effbTag; double weightBTag;double totalWeight;
   for(std::vector<pat::Jet>::const_iterator itJet = myJets->begin(); itJet != myJets->end(); ++itJet){
        newJet = *itJet;
	jetFlav = fabs(itJet->partonFlavour());
	jetPt = itJet->pt() >= 1000.0 ? 997. : itJet->pt() ; jetEta = fabs(itJet->eta()) >= 2.4 ? 2.399 : fabs(itJet->eta());
        if(jetFlav == 5){
          jetPt = jetPt >= 200.0 ? 199.0 : jetPt;
	  effbTag = h_eff_bq->GetBinContent( h_eff_bq->FindBin(jetPt,jetEta) );
	  weightBTag = 1.01;
	  //std::cout<<"h_eff_bq pt "<< jetPt <<" jetEta "<<jetEta<<" eff "<< effbTag <<std::endl;
	}
	else if(jetFlav == 4){
	  jetPt = jetPt >= 200.0 ? 199.0 : jetPt;
	  effbTag = h_eff_cq->GetBinContent( h_eff_cq->FindBin(jetPt,jetEta) );
	  //std::cout<<" h_eff_cq pt "<< jetPt <<" jetEta "<<jetEta<<" eff "<< effbTag <<std::endl;
	  weightBTag = 1.01;
	}
	else {
		//jetPt = jetPt >= 520 ? 519 : jetPt;
		effbTag = h_eff_lq->GetBinContent( h_eff_lq->FindBin(jetPt,jetEta) );
		edm::LogInfo("weightBtag") <<"jetPt"<< jetPt <<"will be reduced to 199";
		weightBTag = btagHist->GetBinContent( btagHist->FindBin(jetPt,jetEta) ) ;// std::cout<<"weightBTag btagHist- "<<weightBTag<<std::endl;
		jetPt = jetPt >= 200.0 ? 199.0 : jetPt;
		effbTag = h_eff_lq->GetBinContent( h_eff_lq->FindBin(jetPt,jetEta) );
	}
	edm::LogInfo("weightBtag")<<"weightBTag "<<weightBTag <<" effbTag "<<effbTag<<" jetPt "<<jetPt<<" jetEta "<<jetEta;
	totalWeight = weightBTag * effbTag; totalWeight = totalWeight > 1.0 ? 1.0 : totalWeight; totalWeight = totalWeight < 0.0 ? 0.0 : totalWeight;
	proba_jetb.push_back(totalWeight); edm::LogInfo("weightBtag")<<"Finally pt "<<jetPt<<" totalWeight "<<totalWeight;
        newJet.addUserFloat("probBJet",totalWeight);newJets.push_back(newJet);
   }//
   //for (unsigned int i = 0; i < newJets.size(); i++){std::cout<<" "<<newJets[i].userFloat("probBJet")<<std::endl;}
   std::auto_ptr<std::vector<pat::Jet> > pOut(new std::vector<pat::Jet>(newJets));
   iEvent.put(pOut);
 

   // is put into the Event
   edm::LogInfo("weightBtag")<<"addBTagWeights bTag event weights ";
   btagWeights.resize(5);btagWeights[0]=1;btagWeights[1]=0;btagWeights[2]=0;btagWeights[3]=0;btagWeights[4]=0;
   // calculating probabilities
           float proba_0jet = 1.;
           float proba_1jet = 0.;
           float proba_2jet = 0.;
           float proba_atleast3jet = 0.;
           for (unsigned int ijet=0; ijet<proba_jetb.size(); ijet++) {
               proba_0jet *= (1-proba_jetb[ijet]);
               float pp1 = 1;
               for (unsigned int kjet=0; kjet<proba_jetb.size(); kjet++) {
                  if (kjet!=ijet) {
                      pp1 *= (1-proba_jetb[kjet]);
                   }

                   if (kjet>ijet) {
                       float pp2 = 1;
                       for (unsigned int ljet=0; ljet<proba_jetb.size(); ljet++) {
                          if (ljet!=ijet && ljet!=kjet) {
                             pp2 *= (1-proba_jetb[ljet]);
                           }
                        }
                        proba_2jet+=proba_jetb[ijet]*proba_jetb[kjet]*pp2;
                   }
               }
               proba_1jet+=proba_jetb[ijet]*pp1;
            }
            if(proba_jetb.size()>2) proba_atleast3jet= 1- proba_0jet - proba_1jet - proba_2jet;
            if (n_bjets_==0) {
              btagWeights[0]*=1.;
            }
            else if (n_bjets_==1) {
              btagWeights[0]*=1 - proba_0jet;
            }
            else if (n_bjets_==2) {
              btagWeights[0]*=1 - proba_0jet - proba_1jet;
            }
            btagWeights[1]=proba_0jet;
            btagWeights[2]=proba_1jet;
            btagWeights[3]=proba_2jet;
            btagWeights[4]=proba_atleast3jet;
	edm::LogInfo("weightBtag")<<"btagWeights[0] "<<btagWeights[0]<<" btagWeights[1] "<<btagWeights[1]<< "  btagWeights[2] "<<btagWeights[2]<<" btagWeights[3] " <<btagWeights[3] << " btagWeights[4] "<<btagWeights[4];

//
//   using namespace edm;
//   std::auto_ptr<std::vector<double> > pOut(new std::vector<double>(btagWeights));
//   iEvent.put(pOut);

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
AddMyBTagWeights::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AddMyBTagWeights::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
AddMyBTagWeights::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
AddMyBTagWeights::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AddMyBTagWeights::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AddMyBTagWeights::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AddMyBTagWeights::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddMyBTagWeights);
