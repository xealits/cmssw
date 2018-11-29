import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.TauDiscriminatorTools import requireLeadTrackCalo

caloRecoTauDiscriminationByECALIsolation = cms.EDProducer("CaloRecoTauDiscriminationByIsolation",

    CaloTauProducer = cms.InputTag('caloRecoTauProducer'),

    Prediscriminants = requireLeadTrackCalo,

    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    TrackerIsolAnnulus_maximumOccupancy   = cms.uint32(0),

    ApplyDiscriminationByECALIsolation    = cms.bool(True),
    ECALisolAnnulus_maximumSumEtCut       = cms.double(1.5)                                                  
)


