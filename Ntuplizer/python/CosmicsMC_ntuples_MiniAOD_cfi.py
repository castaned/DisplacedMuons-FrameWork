import FWCore.ParameterSet.Config as cms

ntuples = cms.EDAnalyzer('ntuplizer',
    nameOfOutput = cms.string('CosmicsMC_Ntuples.root'),
    isData                        = cms.bool(False),
    isAOD                         = cms.bool(False),
    EventInfo                     = cms.InputTag("generator"),
    RunInfo                       = cms.InputTag("generator"),
    BeamSpot                      = cms.InputTag("offlineBeamSpot"),
    displacedGlobalCollection     = cms.InputTag("displacedGlobalMuons"),
    displacedStandAloneCollection = cms.InputTag("displacedStandAloneMuons"),
    displacedMuonCollection       = cms.InputTag("slimmedDisplacedMuons"),
    genParticleCollection         = cms.InputTag("prunedGenParticles"),

    #prescales  = cms.InputTag("patTrigger"),
    bits       = cms.InputTag("TriggerResults","","HLT"),
    #objects    = cms.InputTag("slimmedPatTrigger")
)
