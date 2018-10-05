import FWCore.ParameterSet.Config as cms
#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection ### for updating jets!
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("TagPrepper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#print "isData: ", options.isData
#if options.isData:
    # print "Running on Data"
    # process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#else:
print "Running on MC"
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      'file:../test/myOutputFile_COMB_QG_5ev.root',
    )
)

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)



process.CorrSequence = cms.Sequence(process.patJetCorrFactorsFinal * process.updatedPatJetsFinal)
################################

process.myNtupler = cms.EDProducer('prepper', #src = cms.InputTag("selectedUpdatedPatJetsDeepCSV")
)
process.out = cms.OutputModule("PoolOutputModule",
    #outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('prepOutput.root')
)


process.p = cms.Path(process.BTAGSequence * process.CorrSequence * process.myNtupler)
process.e = cms.EndPath(process.out)
