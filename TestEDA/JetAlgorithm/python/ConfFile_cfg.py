import FWCore.ParameterSet.Config as cms

process = cms.Process("JetFinder")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('L1Trigger.RegionalCaloTrigger.rctDigis_cfi')

process.es_ascii = cms.ESSource("HcalTextCalibrations",
    input = cms.VPSet(
            cms.PSet(
                object = cms.string('LutMetadata'),
                file = cms.FileInPath('CalibCalorimetry/CaloTPG/HcalLutMetadata_now_with_HTV1')
            )
    )
)
process.es_prefer_es_ascii = cms.ESPrefer("HcalTextCalibrations", "es_ascii")

process.boolTrue = cms.EDFilter( 'HLTBool',
    result = cms.bool( True )
)
process.calibPreSequence = cms.Sequence(process.boolTrue)
process.GlobalTag.globaltag = cms.string('GR_R_71_V1::All')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:DIGIoutput.root'
    )
)

process.jetFinder = cms.EDAnalyzer('JetAlgorithm'
)


process.p = cms.Path(process.calibPreSequence + process.jetFinder)
