import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("SLntupler")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
#process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/n/nmangane/LTW3/Demo/prepper/test/BTag_ID_100ev.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1A257ECE-84CF-E611-A64F-141877410512.root',  
    )
)

process.tree = cms.EDAnalyzer('SLntupler')
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('SLntuple.root'),
                                   closeFileFast = cms.untracked.bool(True) 
                                   )
process.p = cms.Path(process.tree)
