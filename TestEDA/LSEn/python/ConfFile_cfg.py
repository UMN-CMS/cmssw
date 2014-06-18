import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.es_ascii = cms.ESSource("HcalTextCalibrations",
    input = cms.VPSet(
        cms.PSet(
            object = cms.string('LutMetadata'),
            file = cms.FileInPath('CalibCalorimetry/CaloTPG/HcalLutMetadata_now_with_HTV1')
        )
    )
)
process.es_prefer_es_ascii = cms.ESPrefer("HcalTextCalibrations", "es_ascii")



process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_71_V1::All')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/E035CA45-6ADC-E111-AF21-BCAEC518FF6B.root','file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/085550A9-3EDC-E111-AA7F-5404A63886BE.root', 'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/08E28448-79DC-E111-9F48-5404A63886C5.root', 'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/18B86D77-5DDC-E111-B6AA-5404A63886EF.root' 
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("testing.root")
)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.demo = cms.EDAnalyzer('LSEn')
process.p = cms.Path(process.demo)
#process.p = cms.Path(process.demo)
process.schedule = cms.Schedule(process.raw2digi_step,process.p)
