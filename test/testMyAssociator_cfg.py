import FWCore.ParameterSet.Config as cms

process = cms.Process("myana")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:RelVal_GEN-SIM-RECO.root'),
    secondaryFileNames = cms.untracked.vstring('file:RelVal_GEN-SIM-DIGI-RAW-HLTDEBUG.root')
)

# MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories = cms.untracked.vstring('testMyAssociator',
    'MuonAssociatorByHits', 'DTHitAssociator', 'RPCHitAssociator', 'MuonTruth', 'PSimHitMap',
    'FwkJob', 'FwkReport', 'FwkSummary', 'Root_NoDictionary')

process.MessageLogger.cout = cms.untracked.PSet(
    noTimeStamps = cms.untracked.bool(True),
    threshold = cms.untracked.string('INFO'),
    INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    testMyAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    ),
    MuonAssociatorByHits = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    DTHitAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    RPCHitAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    MuonTruth = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    PSimHitMap = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    FwkReport = cms.untracked.PSet(
        reportEvery = cms.untracked.int32(1),
        limit = cms.untracked.int32(10000000)
    ),
    FwkSummary = cms.untracked.PSet(
        reportEvery = cms.untracked.int32(1),
        limit = cms.untracked.int32(10000000)
    ),
    FwkJob = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    Root_NoDictionary = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    )
)

process.MessageLogger.statistics = cms.untracked.vstring('cout')

process.MessageLogger.cerr = cms.untracked.PSet(
    placeholder = cms.untracked.bool(True)
)

#process.Tracer = cms.Service("Tracer")

# Mixing Module
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

# Standard Sequences
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_30X::All')

# --- example Analyzer running MuonAssociatorByHits 
process.testanalyzer = cms.EDAnalyzer("testMyAssociator",
    # for Muon Track association
    #
    #     input collections
    #
    # ... reco::Track collection
    tracksTag = cms.InputTag("standAloneMuons"),
    # tracksTag = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    # tracksTag = cms.InputTag("globalMuons"),
    # tracksTag = cms.InputTag("generalTracks"),
    #
    # ... TrackingParticle collection
    tpTag = cms.InputTag("mergedtruth","MergedTrackTruth"),
    #
    dumpInputCollections = cms.bool(False),
    #
    #....... general input parameters
    #
    AbsoluteNumberOfHits_track = cms.bool(False),
    MinHitCut_track = cms.uint32(1),
    AbsoluteNumberOfHits_muon = cms.bool(False),
    MinHitCut_muon = cms.uint32(1),
    #
    UseTracker = cms.bool(False),
    UseMuon = cms.bool(True),
    #
    PurityCut_track = cms.double(0.5),
    PurityCut_muon = cms.double(0.5),
    #
    EfficiencyCut_track = cms.double(0.5),
    EfficiencyCut_muon = cms.double(0.5),
    #
    #........(for inner tracker stub of Global Muons)...
    UsePixels = cms.bool(True),
    UseGrouped = cms.bool(True),
    UseSplitting = cms.bool(True),
    ThreeHitTracksAreSpecial = cms.bool(True),
    #
    # for DT Hit associator
    crossingframe = cms.bool(True),
    simtracksTag = cms.InputTag("g4SimHits"),
    simtracksXFTag = cms.InputTag("mix","g4SimHits"),
    #
    DTsimhitsTag = cms.InputTag("g4SimHits","MuonDTHits"),
    DTsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonDTHits"),
    DTdigiTag = cms.InputTag("simMuonDTDigis"),
    DTdigisimlinkTag = cms.InputTag("simMuonDTDigis"),
    DTrechitTag = cms.InputTag("dt1DRecHits"),
    #
    dumpDT = cms.bool(False),
    links_exist = cms.bool(True),
    associatorByWire = cms.bool(False),
    #
    # for CSC Hit associator
    CSCsimHitsXFTag = cms.InputTag("mix","g4SimHitsMuonCSCHits"),
    CSClinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
    CSCwireLinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigiSimLinks"),
    #
    # for RPC Hit associator
    RPCsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonRPCHits"),
    RPCdigisimlinkTag = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"),
    #
    # for Tracker Hit associator
    #
    associatePixel = cms.bool(True),
    associateStrip = cms.bool(True),
    associateRecoTracks = cms.bool(True),
    #                                
    ROUList = cms.vstring('TrackerHitsTIBLowTof', 
        'TrackerHitsTIBHighTof', 
        'TrackerHitsTIDLowTof', 
        'TrackerHitsTIDHighTof', 
        'TrackerHitsTOBLowTof', 
        'TrackerHitsTOBHighTof', 
        'TrackerHitsTECLowTof', 
        'TrackerHitsTECHighTof', 
        'TrackerHitsPixelBarrelLowTof', 
        'TrackerHitsPixelBarrelHighTof', 
        'TrackerHitsPixelEndcapLowTof', 
        'TrackerHitsPixelEndcapHighTof')
)

process.test = cms.Path(process.mix * process.testanalyzer)
