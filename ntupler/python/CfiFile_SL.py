import FWCore.ParameterSet.Config as cms
tree = cms.EDAnalyzer('SLntupler',
                              isData = cms.bool(False),
                              isMC = cms.bool(False),
                              deBug = cms.bool(False),
                              is2016 = cms.bool(False),
                              is2017 = cms.bool(False),
                              is2018 = cms.bool(False),
                              HTMin = cms.double(0),
                              eleVetoIdMap2016 = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                              eleLooseIdMap2016 = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                              eleMediumIdMap2016 = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                              eleTightIdMap2016 = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
                              )
