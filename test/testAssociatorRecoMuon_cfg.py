import FWCore.ParameterSet.Config as cms

process = cms.Process("myana")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

from Configuration.EventContent.EventContent_cff import *
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'rfio:/castor/cern.ch/user/g/gpetrucc/900GeV/MC/Feb9Skims/MC_CollisionEvents_MuonSkim_1.root'
        'file:/data/gpetrucc/Feb9Skims/MC_CollisionEvents_MuonSkim_1.root'
    ),
    inputCommands = RECOSIMEventContent.outputCommands,            # keep only RECO out of RAW+RECO, for tests
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),  # keep only RECO out of RAW+RECO, for tests
)

# MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories = cms.untracked.vstring('testAssociatorRecoMuon', 'MuonAssociatorByHits')
process.MessageLogger.cout = cms.untracked.PSet(
    noTimeStamps = cms.untracked.bool(True),
    threshold = cms.untracked.string('INFO'),
    INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    default = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    testAssociatorRecoMuon = cms.untracked.PSet(limit = cms.untracked.int32(10000000))
)
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))

# Mixing Module
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

# Standard Sequences
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")           # On RAW+RECO
#process.load("SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi")          # On RAW+RECO
process.load("SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi")   # On RECO
process.load("SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi") # On RECO

process.GlobalTag.globaltag = cms.string('START3X_V26::All')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# --- example Analyzer running MuonAssociatorByHits 
process.testanalyzer = cms.EDAnalyzer("testAssociatorRecoMuon",
    muonsTag  = cms.InputTag("muons"),
    trackType = cms.string("segments"),  # or 'inner','outer','global'
    tpTag    = cms.InputTag("mergedtruthNoSimHits","MergedTrackTruth"),  # RECO Only
    associatorLabel = cms.string("muonAssociatorByHits_NoSimHits"),      # RECO Only
    #tpTag    = cms.InputTag("mergedtruth","MergedTrackTruth"),          # RAW+RECO
    #associatorLabel = cms.string("muonAssociatorByHits"),               # RAW+RECO
) 

process.skim = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("muons"), minNumber = cms.uint32(1))
process.test  = cms.Path(process.skim+process.mix * process.trackingParticlesNoSimHits * process.testanalyzer) # RECO
#process.test = cms.Path(process.skim+process.mix * process.trackingParticles          * process.testanalyzer) # RAW+RECO
