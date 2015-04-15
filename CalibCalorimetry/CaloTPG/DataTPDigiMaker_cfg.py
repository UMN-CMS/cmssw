import FWCore.ParameterSet.Config as cms 
  
process = cms.Process("DataDIGIs")


#--- Steps available to this configuration       ---#
#--- (1) Trigger filter.  Input data RAW/AOD     ---#
#--- (2) N(vertex) filter.  Only valid for MC    ---#
#--- (3) Reconstruction, assuming RAW input data ---#
#--- (4) Filtering on ECAL+HF, Z mass            ---#
#--- (5) HF Calibration analysis                 ---#
doRerecoOnRaw = True

#-----------------------------------#
#--- Flag for running on data/MC ---#
#-----------------------------------#
isData = True

## Import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

process.load('L1Trigger.Configuration.L1Extra_cff')

#############added

#################3added

## Global tags:
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# Note that this assumes the 6X post-LS1 Monte Carlo

if( isData):
    process.GlobalTag.globaltag = cms.string('GR_R_71_V1::All')
#####################fix this
else:
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')


process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigiData_cff')
#process.HcalTrigTowerGeometryESProducer.useFullGranularityHF = cms.bool( True )



process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

if(isData):
      TempFN=['file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/085550A9-3EDC-E111-AA7F-5404A63886BE.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/08E28448-79DC-E111-9F48-5404A63886C5.root', 
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/18B86D77-5DDC-E111-B6AA-5404A63886EF.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/00484FFF-4EDC-E111-8A8F-0025901D629C.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/1AB266F2-56DC-E111-BBCD-0025901D631E.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/1AE9FEE6-2EDC-E111-82E9-E0CB4E55365D.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/1EFD93DA-33DC-E111-AADC-BCAEC532971B.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/20AC7ED2-3BDC-E111-BD44-BCAEC518FF8A.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/20BA83FE-38DC-E111-A63A-5404A63886C5.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/30FE52F7-44DC-E111-9487-BCAEC518FF50.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/326EA6E8-4CDC-E111-8C2F-5404A63886B9.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/00484FFF-4EDC-E111-8A8F-0025901D629C.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/340603FF-49DC-E111-8283-BCAEC5364C62.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/36A64EDB-5CDC-E111-99C7-001D09F29619.root',
      'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/40747533-33DC-E111-A167-003048D373AE.root'
      ]
      #TempFN=['file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/201/278/02E206EC-31EB-E111-B28B-0025901D625A.root']
              
      
else:
    TempFN=['file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/48CC1D33-6E6A-E311-8855-001D09F251EF.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/5668ACA2-6B6A-E311-A71E-02163E00E76C.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/5E0D6BBD-6C6A-E311-9202-0025904B2FD8.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/608D9D23-6B6A-E311-AEC2-E0CB4E55365D.root']


process.source = cms.Source("PoolSource",
    # eventsToProcess = cms.untracked.VEventRange(),
    #fileNames = cms.untracked.vstring( 'file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/48CC1D33-6E6A-E311-8855-001D09F251EF.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/5668ACA2-6B6A-E311-A71E-02163E00E76C.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/5E0D6BBD-6C6A-E311-9202-0025904B2FD8.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/608D9D23-6B6A-E311-AEC2-E0CB4E55365D.root')
    
fileNames = cms.untracked.vstring(TempFN )
    
)

process.out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string("/data/whybee0b/user/aevans/DataDIGIoutput.root"),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring( 'keep *_*simHcalTriggerPr*_*_*','keep *_*_*Forward*_*')
    
)

###--- (3) Re-RECO from RAW ---###
### Auto generated configuration file using Revision: 1.381.2.6 
### Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
### with command line options: -s RAW2DIGI,RECO ...

if doRerecoOnRaw:
    # process.reconstructionFromRawSequence = cms.Sequence(process.RawToDigi * process.L1Reco * process.reconstruction)
    process.reconstructionFromRawSequence = cms.Sequence(process.RawToDigi)

# From Maria Cepeda Hermida's suggestion for rctDigis testing
#from L1Trigger.RegionalCaloTrigger.rctDigis_cfi import rctDigis
#process.rctDigis = rctDigis
process.load('L1Trigger.RegionalCaloTrigger.rctDigis_cfi')
process.rctDigis.hcalDigis = cms.VInputTag(cms.InputTag("simHcalTriggerPrimitiveDigis"))


#--- Create TP digis from unpacked digis
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cfi import LSParameter
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalUnsuppressedDigis'), cms.InputTag('simHcalUnsuppressedDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

#--- Conditions dump ---#
#process.hcalCond = cms.EDAnalyzer("HcalDumpConditions",
#    dump = cms.untracked.vstring(
#         'LutMetadata'
#    ),
#    outFilePrefix = cms.untracked.string('DumpCond')
#)

process.es_ascii = cms.ESSource("HcalTextCalibrations",
    input = cms.VPSet(
        cms.PSet(
            object = cms.string('LutMetadata'),
            file = cms.FileInPath('CalibCalorimetry/CaloTPG/HcalLutMetadata_now_with_HTV1')
        )
    )
)

process.lumiProducer = cms.EDProducer('LumiProducer',
  connect = cms.string('frontier://LumiProd/CMS_LUMI_PROD'),
  lumiversion = cms.untracked.string(''),
  ncacheEntries = cms.untracked.uint32(5),
)

LumiCorrectionSource = cms.ESSource( "LumiCorrectionSource",
  connect = cms.string('frontier://LumiCalc/CMS_LUMI_PROD'),
)
process.es_prefer_es_ascii = cms.ESPrefer("HcalTextCalibrations", "es_ascii")

#####3333
#process.source= cms.Source("PoolSource",
#              processingMode=cms.untracked.string('RunsAndLumis')        
              #fileNames=cms.untracked.vstring('file:/data/cmsdata/200786/FC7E661B-C3E8-E111-A23E-003048D2BDD8.root') 
#        )
process.LumiCorrectionSource=cms.ESSource("LumiCorrectionSource",
              connect=cms.string('frontier://LumiCalc/CMS_LUMI_PROD'),
              normtag=cms.untracked.string('HFV2a')    
        )


#--- Dump digi information ---#
process.dump  = cms.EDAnalyzer("HcalDigiDump")
process.test = cms.EDAnalyzer("tester")
process.DataAna =cms.EDAnalyzer("DataAnalizer")
process.CutAs = cms.EDAnalyzer("CutAssessment")
process.FullGraph = cms.EDAnalyzer('LSEn')
#--- Histograms ---#
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cfi import LSParameter
#new_filename = "Long"+LSParameter.Min_Long_Energy.pythonValue() + "_Short" + LSParameter.Min_Short_Energy.pythonValue() + "_Slope_OffSet" + LSParameter.Long_vrs_Short_Slope.pythonValue() + "_Long_Short_Offset"+LSParameter.Long_Short_Offset.pythonValue() +".root"
new_filename = "Long%s_Short%s_Slope%s_LongShortOffset%s.root" % (LSParameter.Min_Long_Energy.pythonValue(),
                                                                  LSParameter.Min_Short_Energy.pythonValue(),
                                                                  LSParameter.Long_vrs_Short_Slope.pythonValue(),
                                                                  LSParameter.Long_Short_Offset.pythonValue())
                                                                  
if(isData):
    new_filename="fdata/"+new_filename
else:
    new_filename="fmc/"+new_filename
    new_filename="Datatest.root"
    process.TFileService = cms.Service("TFileService",
    fileName = cms.string(new_filename)
    )
#process.histos = cms.EDAnalyzer("HcalTriggerDigiDump")

###--- Assemble everything ---###
process.boolTrue = cms.EDFilter( 'HLTBool',
    result = cms.bool( True )
)
process.raw2digi_step = cms.Path(process.RawToDigi)

process.calibPreSequence = cms.Sequence(process.boolTrue)
######################
# Additional output definition


# Path and EndPath definitions
#process.digitisation_step = cms.Path(process.pdigi)
#process.L1simulation_step = cms.Path(process.SimL1Emulator)
#process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.out)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
# Schedule definition
#process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step)


###########################
# process.p = cms.Path( process.calibPreSequence + process.hcalCond + process.simHcalTriggerPrimitiveDigis + process.dump + process.histos )
#process.p = cms.Path( process.calibPreSequence + process.simHcalTriggerPrimitiveDigis + process.test + process.CutAs + process.dump + process.rctDigis )

if( isData):
    process.p = cms.Path(process.L1Extra + process.calibPreSequence  + process.simHcalTriggerPrimitiveDigis+process.lumiProducer)
    #process.reconstructionFromRawSequence+
    #process.p = cms.Path(process.reconstructionFromRawSequence+ process.calibPreSequence  + process.simHcalTriggerPrimitiveDigis + process.DataAna+process.FullGraph)
else:
    process.p = cms.Path( process.calibPreSequence + process.simHcalTriggerPrimitiveDigis + process.test + process.CutAs)
#process.RECOSIMoutput_step = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.raw2digi_step,process.p,process.endjob_step)



