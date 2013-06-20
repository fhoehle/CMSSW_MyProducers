// -*- C++ -*-
//
// Package:    AddPileUpWeightsProducer
// Class:      AddPileUpWeightsProducer
// 
/**\class AddPileUpWeightsProducer AddPileUpWeightsProducer.cc CMSSW_MyProducers/AddPileUpWeightsProducer/src/AddPileUpWeightsProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Felix Hoehle
//         Created:  Thu Mar 29 00:39:36 CEST 2012
// $Id: AddPileUpWeightsProducer.cc,v 1.2 2012/08/11 16:34:50 fhohle Exp $
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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
// class declaration
//

class AddPileUpWeightsProducer : public edm::EDProducer {
   public:
      explicit AddPileUpWeightsProducer(const edm::ParameterSet&);
      ~AddPileUpWeightsProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag vertexSrc_;

  edm::LumiReWeighting LumiWeights_;
  std::string PileupFile1_;
  std::string PileupFile2_;
  std::string PUHistname2_;
  std::string PUHistname1_;

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
AddPileUpWeightsProducer::AddPileUpWeightsProducer(const edm::ParameterSet& iConfig):
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
PileupFile1_(iConfig.getParameter<std::string>("pileupFile1")),
PileupFile2_(iConfig.getParameter<std::string>("pileupFile2")),
PUHistname1_(iConfig.getParameter<std::string>("PUHistname1")),
PUHistname2_(iConfig.getParameter<std::string>("PUHistname2"))
{
   //register your products
//   produces<std::vector<std::pair<std::basic_string<char>,std::vector<float> > > >();
 produces<std::vector< double > >();
/*
   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  LumiWeights_ = edm::LumiReWeighting(PileupFile1_.c_str(),PileupFile2_.c_str(),PUHistname1_.c_str(),PUHistname2_.c_str());
  edm::LogInfo   ("lumiweights") << "using files" << "file1" << PileupFile1_.c_str() << "hist file1"<< PUHistname1_.c_str() <<"file2"<< PileupFile2_.c_str() << "hist file2"<< PUHistname2_.c_str();
}


AddPileUpWeightsProducer::~AddPileUpWeightsProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
AddPileUpWeightsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

  std::vector<PileupSummaryInfo>::const_iterator PVI;


  float npT=-1.;
  float npIT=-1.;
  float npBCm1IT = -1.0;
  float npBCp1IT = -1.0;
  //
  ////std::cout << "pu summary 189 " << PupInfo->size() << std::endl;
  //// (then, for example, you can do)
  //for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
  //  int n_bc=PVI->getBunchCrossing();
  //  //std::cout << " Pileup Information: bunchXing, nvtx, size : " << PVI->getBunchCrossing() << " " 
  //  //<< PVI->getPU_NumInteractions()  << " " 
  //  //<< PVI->getPU_zpositions().size() << std::endl;
  //
  //  if (n_bc==0 )  cev.num_pileup_bc0    = PVI->getPU_NumInteractions();
  //  if (n_bc==1 )  cev.num_pileup_bcp1   = PVI->getPU_NumInteractions();
  //  if (n_bc==-1)  cev.num_pileup_bcm1   = PVI->getPU_NumInteractions();
  //
  //
  //}

  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

    int BX = PVI->getBunchCrossing();
    if (BX == -1) npBCm1IT = PVI->getPU_NumInteractions();
    if(BX == 0) {
      npT = PVI->getTrueNumInteractions();
      npIT = PVI->getPU_NumInteractions();
    }
    if (BX == 1) npBCp1IT = PVI->getPU_NumInteractions();
  }
  edm::LogInfo   ("lumiweights") << "before calling LumiWeights";
  double MyWeight_npT = LumiWeights_.weight( npT );
  double MyWeight_npIT = LumiWeights_.weight( npIT );
  double MyWeightIT_npIT = LumiWeights_.weight(npIT );
  double MyWeight_npavIT = LumiWeights_.weight( float((npBCm1IT+npIT+npBCp1IT)/3.0) );
  edm::LogInfo   ("lumiweights") << "calling LumiWeights done";
//   using namespace edm;
   std::vector<double> weights; 
   edm::LogInfo   ("lumiweights") <<iEvent.id ().run()<<"  lumi "<< iEvent.id ().luminosityBlock ()<<" event "<<iEvent.id ().event()<<" MyWeight_npT "<<MyWeight_npT<<std::endl; weights.push_back(MyWeight_npT);
   edm::LogInfo   ("lumiweights")<<iEvent.id ().run()<<"  lumi "<< iEvent.id ().luminosityBlock ()<<" event "<<iEvent.id ().event()<<" MyWeight_npIT "<<MyWeight_npIT<<std::endl;  weights.push_back(MyWeight_npIT);
   edm::LogInfo   ("lumiweights")<<iEvent.id ().run()<<"  lumi "<< iEvent.id ().luminosityBlock ()<<" event "<<iEvent.id ().event()<<" MyWeightIT_npIT "<<MyWeightIT_npIT<<std::endl; weights.push_back(MyWeightIT_npIT);
   edm::LogInfo   ("lumiweights")<<"PU weights "; for(unsigned int i = 0;i < weights.size();i++){edm::LogInfo   ("lumiweights")<<" "<<weights[i];} edm::LogInfo   ("lumiweights")<<" end "<<std::endl;
   std::auto_ptr< std::vector< double > > pOut (new std::vector<double>(weights));
   edm::LogInfo   ("lumiweights")<<"PU pOut weights "; for(unsigned int i = 0;i < pOut->size();i++){std::cout<<" "<<(*pOut)[i];} std::cout<<" end "<<std::endl;
   iEvent.put(pOut);


/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
AddPileUpWeightsProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AddPileUpWeightsProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
AddPileUpWeightsProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
AddPileUpWeightsProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AddPileUpWeightsProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AddPileUpWeightsProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AddPileUpWeightsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddPileUpWeightsProducer);
