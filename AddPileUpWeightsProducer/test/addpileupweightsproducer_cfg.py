import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:myfile.root'
    )
)

process.myProducerLabel = cms.EDProducer("AddPileUpWeightsProducer",
                                  vertexSrc = cms.InputTag("offlinePrimaryVertices"),
  pileupFile1 = cms.string("$CMSSW_BASE/src/CMSSW_CMSSW_MyProducers/AddPileUpWeightsProducer/data/pileUpHistos_TTbar-SemiMu_Selected_Reconstructed.root"),
  pileupFile2 = cms.string("$CMSSW_BASE/src/CMSSW_CMSSW_MyProducers/AddPileUpWeightsProducer/data/Test_Data_Run2011A_Pileup.root"),
  PUHistname1 = cms.string("h1_TNPUTrue"),
  PUHistname2 = cms.string("pileup") 
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
