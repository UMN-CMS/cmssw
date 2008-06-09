# The following comments couldn't be translated into the new config version:

# All/OuterSurface/InnerSurface/ImpactPoint/default(track)
#

import FWCore.ParameterSet.Config as cms

TrackMon = cms.EDFilter("TrackingMonitor",
    OutputMEsInRootFile = cms.bool(False),
    phiErrMax = cms.double(5.0),
    MeasurementState = cms.string('default'),
    etaErrMax = cms.double(0.5),
    TrackPxBin = cms.int32(100),
    RecHitBin = cms.int32(22),
    TrackPzMin = cms.double(-50.0),
    Chi2Max = cms.double(500.0),
    Chi2Bin = cms.int32(100),
    TrackPzBin = cms.int32(100),
    pxErrBin = cms.int32(1000),
    etaErrMin = cms.double(-0.5),
    TrackPxMax = cms.double(50.0),
    TrackPzMax = cms.double(50.0),
    ThetaBin = cms.int32(100),
    RecHitMin = cms.double(0.0),
    OutputFileName = cms.string('MonitorTrack.root'),
    pzErrMin = cms.double(-1000.0),
    EtaMin = cms.double(-4.0),
    pErrBin = cms.int32(1000),
    pyErrMin = cms.double(-1000.0),
    phiErrBin = cms.int32(500),
    EtaMax = cms.double(4.0),
    etaErrBin = cms.int32(500),
    Chi2Min = cms.double(-0.5),
    ThetaMin = cms.double(0.0),
    PhiMin = cms.double(-3.2),
    TrackPtMax = cms.double(30.0),
    RecHitMax = cms.double(25.0),
    TrackPyBin = cms.int32(100),
    ptErrMin = cms.double(0.0),
    TkSizeMin = cms.double(0.0),
    TrackPxMin = cms.double(-50.0),
    pyErrMax = cms.double(1000.0),
    ThetaMax = cms.double(3.2),
    pzErrMax = cms.double(1000.0),
    pzErrBin = cms.int32(1000),
    pErrMin = cms.double(-1000.0),
    EtaBin = cms.int32(100),
    pErrMax = cms.double(1000.0),
    #
    FolderName = cms.string('Track/GlobalParameters'),
    pxErrMax = cms.double(1000.0),
    TkSizeBin = cms.int32(500),
    AlgoName = cms.string('GenTk'),
    TrackPyMin = cms.double(-50.0),
    TrackProducer = cms.InputTag("generalTracks"),
    TkSizeMax = cms.double(500.0),
    TrackPtBin = cms.int32(100),
    TrackPyMax = cms.double(50.0),
    phiErrMin = cms.double(-5.0),
    pyErrBin = cms.int32(1000),
    pxErrMin = cms.double(-1000.0),
    PhiBin = cms.int32(100),
    ptErrMax = cms.double(1000.0),
    PhiMax = cms.double(3.2),
    TrackPtMin = cms.double(-0.5),
    ptErrBin = cms.int32(500)
)



