import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/Fall11/TT_TuneZ2_7TeV-mcatnlo/AODSIM/PU_S6_START42_V14B-v1/0000/0285DC46-342A-E111-BCE0-00304867D838.root'
    )
)
process.maxEvents.input=10
process.myProducerLabel = cms.EDProducer('MyTTbarGenEvent10Parts'
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printList = cms.EDAnalyzer("ParticleListDrawer",
        src = cms.InputTag("genParticles"),
        maxEventsToPrint = cms.untracked.int32(10)
)

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
        src = cms.InputTag("genParticles"),
        printP4 = cms.untracked.bool(True),
        printPtEtaPhi = cms.untracked.bool(False),
        printVertex = cms.untracked.bool(False),#True),
        printStatus = cms.untracked.bool(False),
        printIndex = cms.untracked.bool(False),
        status = cms.untracked.vint32(3)#, 2, 3)
)
  
process.p = cms.Path(process.myProducerLabel * process.printList )

process.e = cms.EndPath(process.out)
