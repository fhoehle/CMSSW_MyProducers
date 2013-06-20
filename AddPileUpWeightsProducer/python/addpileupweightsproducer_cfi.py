import FWCore.ParameterSet.Config as cms

addPileupWeightsProducer = cms.EDProducer("AddPileUpWeightsProducer",
                                  vertexSrc = cms.InputTag("offlinePrimaryVertices"),
  pileupFile1 = cms.string("$CMSSW_BASE/src/CMSSW_CMSSW_MyProducers/AddPileUpWeightsProducer/input/pileUpHistos_TTbar-SemiMu_Selected_Reconstructed.root"),
  pileupFile2 = cms.string("$CMSSW_BASE/src/CMSSW_CMSSW_MyProducers/AddPileUpWeightsProducer/input/Test_Data_Run2011A_Pileup.root"),
  PUHistname1 = cms.string("h1_TNPUTrue"),
  PUHistname2 = cms.string("pileup") 
)


