import FWCore.ParameterSet.Config as cms
MyTTbarGenEventProd = cms.EDProducer('MyTTbarGenEventLepHad',
        genTag = cms.untracked.InputTag("genParticles"),
        )

