import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )                    

process.load('FWCore.MessageService.MessageLogger_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))          
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(50),SkipEvent = cms.untracked.vstring('ProductNotFound'))


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
        #'/store/data/Run2017C/Charmonium/MINIAOD/PromptReco-v1/000/299/369/00000/18BAE1A7-936D-E711-A321-02163E01A5E3.root',
        '/store/data/Run2017B/Charmonium/MINIAOD/PromptReco-v1/000/297/050/00000/183B4680-3356-E711-B33A-02163E014487.root',

    )
)

process.demo = cms.EDAnalyzer('DemoTrackAnalyzer',
               tracks = cms.InputTag("packedPFCandidates"),
)

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('outfile.root')
                                   )

process.p = cms.Path(process.demo)
