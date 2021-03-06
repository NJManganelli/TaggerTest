import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection ### for updating jets!
from Configuration.AlCa.GlobalTag import GlobalTag

#process = cms.Process("BTagIDProcess")
process = cms.Process("SLntupler")

#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

################################
################################
##### Input File Selection #####
################################
################################
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/FE78DF14-9A33-E811-8CC3-0025907D2000.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/FCF4D5A8-7530-E811-AD27-FA163E34AA7F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/FAA5C1A7-D130-E811-BBC8-7845C4FC379D.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/F091C785-9A33-E811-A2A3-FA163E726F24.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/EE1A4FEE-5132-E811-A2E4-002590E7E012.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/E8E62005-AF30-E811-B4E4-02163E012E9A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/E8A15EA1-9933-E811-925F-008CFAF74A32.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/E891CED7-9933-E811-A1DC-001F2908EC96.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/E4EE72D1-AA30-E811-87B8-FA163E0D52BF.root',
# 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/E07F5AA7-9933-E811-B1B8-0CC47A5FBE31.root',
# 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/DE44CE94-9F31-E811-B7E4-A4BF01025C02.root',
# 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/DCC19830-3B32-E811-B51C-0242AC130002.root',
# 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/DCA1CD0F-1A32-E811-A82C-A4BF011595FC.root',
# 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/D47374EF-5131-E811-A4B0-D4AE52E7FBCA.root',
# 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTToSemiLepton_HT500Njet9_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/D030ACD6-7832-E811-A4E3-1CB72C0A3D89.root',
    )
)

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

################################################
################################################
##### JET Energy Corrections and B-Tagging #####
################################################
################################################

# Configure QGTagger and prepare to add DeepCSV variables #
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag("slimmedJets")
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
process.BTAGSequence = cms.Sequence(process.QGTagger *
                                    process.patJetCorrFactorsDeepCSV * 
                                    process.updatedPatJetsDeepCSV * 
                                    process.pfImpactParameterTagInfosDeepCSV *
                                    process.pfInclusiveSecondaryVertexFinderTagInfosDeepCSV *
                                    process.pfDeepCSVTagInfosDeepCSV *
                                    process.pfDeepCSVJetTagsDeepCSV *
                                    process.patJetCorrFactorsTransientCorrectedDeepCSV *
                                    process.updatedPatJetsTransientCorrectedDeepCSV * 
                                    process.QGTagger *
                                    process.selectedUpdatedPatJetsDeepCSV)


#####################################
#####################################
##### Electron ID Configuration #####
#####################################
#####################################
#Mirena's code

#Electron ID configuration
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff', # previously used with 2016 data
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff', # previously used with 2016 data
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff', # previously used with 2016 data
    # 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff', # with 2017 data
    # 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', # with 2017 data 
    # 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff', # with 2017 data
    ]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)    

# Example Ele Id process configuration #
process.eleIDntupler = cms.EDAnalyzer(
    'ElectronNtuplerVIDDemo',
    # The module automatically detects AOD vs miniAOD, so we configure both
    #
    # Common to all formats objects
    #
###    beamSpot = cms.InputTag('offlineBeamSpot'),
    #
    # Objects specific to AOD format
    #
###    electrons    = cms.InputTag("gedGsfElectrons"),
###    genParticles = cms.InputTag("genParticles"),
###    vertices     = cms.InputTag("offlinePrimaryVertices"),
###    conversions  = cms.InputTag('allConversions'),
    #
    # Objects specific to MiniAOD format
    #
    electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
    genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
    verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),
    #
    # ID decisions (common to all formats)
    #
    eleIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"),
    # An example of configuration for accessing the full cut flow info for
    # one case is shown below.
    # The map name for the full info is the same as the map name of the
    # corresponding simple pass/fail map above, they are distinguished by
    # the type of the content.
    eleIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"),
    # This is a fairly verbose mode if switched on, with full cut flow 
    # diagnostics for each candidate. Use it in a low event count test job.
    eleIdVerbose = cms.bool(False)
    )

#################################
#################################
##### NTupler Configuration #####
#################################
#################################

# Placeholder for integrating ../../ntupler/plugins/ntupler.cc #
process.load("Demo.ntupler.CfiFile_SL")
process.tree.theProblemEvent = cms.uint32(0)
process.tree.isMC = cms.bool(True)
process.tree.is2016 = cms.bool(True)
process.tree.deBug = cms.bool(False)
process.tree.verBose = cms.bool(False)
process.tree.maskDeepCSV = cms.bool(False)
process.tree.HTMin = cms.double(500)
process.tree.NjMin = cms.int32(7)
process.tree.NbMin = cms.int32(2)
# process.myNtupler = cms.EDProducer('prepper', 
#                                    src = cms.InputTag("selectedUpdatedPatJetsDeepCSV"),
# )

#######################
#######################
##### Output File #####
#######################
#######################
#for outputting a full MiniAOD file
# process.out = cms.OutputModule("PoolOutputModule",
#     #outputCommands = cms.untracked.vstring('keep *'),
#     #fileName = cms.untracked.string('myOutputFile_COMB_QGadded_100ev.root')
#     fileName = cms.untracked.string('TEMPNAME.root')
# )

#for the ntupler
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('LongTTSLFilt500HT_corr.root'),
                                   closeFileFast = cms.untracked.bool(True) 
                                   )

#process.p = cms.Path(process.BTAGSequence * process.CorrSequence) # * process.myNtupler)
process.p = cms.Path(process.egmGsfElectronIDSequence * process.BTAGSequence * process.tree)
#process.e = cms.EndPath(process.out) #not needed since ntupler writes directly into file
