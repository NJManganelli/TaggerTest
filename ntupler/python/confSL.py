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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/n/nmangane/LTW3/Demo/prepper/test/BTag_ID_100ev.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1A257ECE-84CF-E611-A64F-141877410512.root',  
    )
)

#Create process.tree with default values. 
#One of isData and isMC must be set to true.
#One of is2016, is2017, is2018 must be set to true.
#deBug may be set to true to dump event details
#HTMin may be set to any double value greater than 0
#NjMin may be set to an int32 value noting minimum number of jets. Less than 0 indicates no cuts
process.load("Demo.ntupler.CfiFile_SL")
process.tree.theProblemEvent = cms.uint32(0)
process.tree.isMC = cms.bool(True)
process.tree.is2016 = cms.bool(True)
process.tree.deBug = cms.bool(False)
process.tree.verBose = cms.bool(False)
process.tree.maskDeepCSV = cms.bool(True)
process.tree.HTMin = cms.double(50)
process.tree.NjMin = cms.int32(3)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('SLntupleMasked.root'),
                                   closeFileFast = cms.untracked.bool(True) 
                                   )
process.p = cms.Path(process.tree)
