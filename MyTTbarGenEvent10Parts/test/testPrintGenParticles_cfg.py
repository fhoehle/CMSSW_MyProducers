import FWCore.ParameterSet.Config as cms

process = cms.Process("testPrintGenParticles")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Fall11/TT_TuneZ2_7TeV-mcatnlo/AODSIM/PU_S6_START42_V14B-v1/0000/98CBA318-412A-E111-ABB8-002618943948.root')
)

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(100)
)


process.printEventNumber = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.printTree)
process.outpath = cms.EndPath(process.printEventNumber)
process.MessageLogger.destinations = cms.untracked.vstring('cout','cerr')
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('eventsToProcess',
                   '',
                   VarParsing.multiplicity.list,
                   VarParsing.varType.string,
                   "Events to process")
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,'skipEvents')
process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
options.parseArguments()
print "default source ",process.source.fileNames
if options.inputFiles != cms.untracked.vstring():
 process.source.fileNames=options.inputFiles
if options.eventsToProcess:
 process.source.eventsToProcess = cms.untracked.VEventRange (options.eventsToProcess)
if options.maxEvents != '':
 print "maxEvents",options.maxEvents
 process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))



