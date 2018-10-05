import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection ### for updating jets!
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("TagPrepper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 20
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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       # 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root',
      'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1A257ECE-84CF-E611-A64F-141877410512.root',  
    )
)

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)
#########################
jetTag = cms.InputTag("slimmedJets")
#First run the QGTagger, so that its output may be added to the jet collection
#Correct jets for the QGTagger, prior to running DeepCSV which uncorrects the jets
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag("slimmedJets")
#process.QGTagger.jec
updateJetCollection(
   process,
   labelName = "DeepCSV",
   jetSource = cms.InputTag("slimmedJets"),
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
   btagDiscriminators = [
      'pfDeepCSVJetTags:probudsg',
      'pfDeepCSVJetTags:probb',
      'pfDeepCSVJetTags:probc',
      'pfDeepCSVJetTags:probbb',
      'pfDeepCSVJetTags:probcc',
      ] ## to add discriminators                                                                                                                                                                                                              
)
process.updatedPatJetsDeepCSV.userData.userFloats.src += ['QGTagger:qgLikelihood', 'QGTagger:ptD', 'QGTagger:axis1', 'QGTagger:axis2',]
process.updatedPatJetsDeepCSV.userData.userInts.src += ['QGTagger:mult',]
process.jecSequence = cms.Sequence(process.patJetCorrFactorsDeepCSV * process.QGTagger * process.updatedPatJetsDeepCSV * process.patJetCorrFactorsDeepCSV)
#process.jecSequence = cms.Sequence(process.patJetCorrFactorsDeepCSV * process.updatedPatJetsDeepCSV)
################################

process.myProducerLabel = cms.EDProducer('prepper', #src = cms.InputTag("selectedUpdatedPatJetsDeepCSV")
)
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('myOutputFile_wTI_reCorr.root')
)


process.p = cms.Path(process.jecSequence * process.myProducerLabel)
process.e = cms.EndPath(process.out)
