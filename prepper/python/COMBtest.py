import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection ### for updating jets!
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
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root',
      'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1A257ECE-84CF-E611-A64F-141877410512.root',  
      #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ttZJets_13TeV_madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/00454275-4BFC-E611-A3F9-02163E01A211.root',
      #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/F2E6D27E-FDFD-E611-AF19-00266CF9B878.root',
    )
)

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

#########################
jetTag = cms.InputTag("slimmedJets")
########## Configure QGTagger and prepare to add DeepCSV variables ###################

#process.load('RecoJets.JetProducers.QGTagger_cfi')
#process.QGTagger.srcJets   = cms.InputTag("slimmedJets")
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
#process.updatedPatJetsDeepCSV.userData.userFloats.src += ['QGTagger:qgLikelihood', 'QGTagger:ptD', 'QGTagger:axis1', 'QGTagger:axis2',]
#process.updatedPatJetsDeepCSV.userData.userInts.src += ['QGTagger:mult',]
#process.BTAGSequence = cms.Sequence(process.patJetCorrFactorsDeepCSV * process.updatedPatJetsDeepCSV * patJetCorrFactorsTransientCorrectedDeepCSV * process.updatedPatJetsTransientCorrectedDeepCSV) #process.QGTagger *
process.BTAGSequence = cms.Sequence(process.patJetCorrFactorsDeepCSV * process.updatedPatJetsDeepCSV * process.selectedUpdatedPatJetsDeepCSV) #process.QGTagger *

########## Recorrect Jets since DeepCSV breaks this in CMSSW_80X##############
#process.load('RecoJets.JetProducers.QGTagger_cfi')
#process.QGTagger.srcJets   = cms.InputTag("updatedPatJetsFinal")
updateJetCollection(
   process,
   labelName = "Final",
   jetSource = cms.InputTag("updatedPatJetsDeepCSV"),
   #jetSource = cms.InputTag("selectedUpdatedPatJetsDeepCSV"), #doesn't exist
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),

)
process.CorrSequence = cms.Sequence(process.patJetCorrFactorsFinal * process.updatedPatJetsFinal)
################################

process.myNtupler = cms.EDProducer('prepper', #src = cms.InputTag("selectedUpdatedPatJetsDeepCSV")
)
process.out = cms.OutputModule("PoolOutputModule",
    #outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('myOutputFile_COMB_QG_100ev.root')
)


#process.p = cms.Path(process.BTAGSequence * process.CorrSequence) # * process.myNtupler)
process.p = cms.Path(process.BTAGSequence)
process.e = cms.EndPath(process.out)
