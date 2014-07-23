import FWCore.ParameterSet.Config as cms
 
process = cms.Process("HFCALIB")


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


process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTrigTowerGeometryESProducer.useFullGranularityHF = cms.bool( True )



process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1000) )

if(isData):
      #TempFN=['file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/085550A9-3EDC-E111-AA7F-5404A63886BE.root','file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/E035CA45-6ADC-E111-AF21-BCAEC518FF6B.root', 'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/08E28448-79DC-E111-9F48-5404A63886C5.root', 'file:/hdfs/cms/phedex/store/data/Run2012C/MinimumBias/RAW/v1/000/200/091/18B86D77-5DDC-E111-B6AA-5404A63886EF.root']
      TempFN=['file:/hdfs/cms/user/lesko/MinBiasWithDigiAndReco/MinBias2012CDigiAndReco_041.root']
      
else:
    TempFN=['file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/48CC1D33-6E6A-E311-8855-001D09F251EF.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/5668ACA2-6B6A-E311-A71E-02163E00E76C.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/5E0D6BBD-6C6A-E311-9202-0025904B2FD8.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/608D9D23-6B6A-E311-AEC2-E0CB4E55365D.root']



process.source = cms.Source("PoolSource",
    # eventsToProcess = cms.untracked.VEventRange(),
    #fileNames = cms.untracked.vstring( 'file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/48CC1D33-6E6A-E311-8855-001D09F251EF.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/5668ACA2-6B6A-E311-A71E-02163E00E76C.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/5E0D6BBD-6C6A-E311-9202-0025904B2FD8.root','file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0_pre11/RelValZEE_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_POSTLS162_V5_OldTrk-v1/00000/608D9D23-6B6A-E311-AEC2-E0CB4E55365D.root')
    
fileNames = cms.untracked.vstring(TempFN )
    
)

#process.out = cms.OutputModule( "PoolOutputModule",
#    fileName = cms.untracked.string("lumprod.root"),
#    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
#    outputCommands = cms.untracked.vstring( 'keep *' )
#)

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
process.calibPreSequence = cms.Sequence(process.boolTrue)
######################
# Additional output definition


# Path and EndPath definitions
#process.digitisation_step = cms.Path(process.pdigi)
#process.L1simulation_step = cms.Path(process.SimL1Emulator)
#process.digi2raw_step = cms.Path(process.DigiToRaw)
#process.endjob_step = cms.EndPath(process.out)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
#process.raw2digi_step = cms.Path(process.RawToDigi)
# Schedule definition
#process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step)


###########################
# process.p = cms.Path( process.calibPreSequence + process.hcalCond + process.simHcalTriggerPrimitiveDigis + process.dump + process.histos )
#process.p = cms.Path( process.calibPreSequence + process.simHcalTriggerPrimitiveDigis + process.test + process.CutAs + process.dump + process.rctDigis )

if( isData):
	process.p = cms.Path(process.calibPreSequence  + process.simHcalTriggerPrimitiveDigis+process.lumiProducer+ process.DataAna)
	#process.reconstructionFromRawSequence+
    #process.p = cms.Path(process.reconstructionFromRawSequence+ process.calibPreSequence  + process.simHcalTriggerPrimitiveDigis + process.DataAna+process.FullGraph)
else:
    process.p = cms.Path( process.calibPreSequence + process.simHcalTriggerPrimitiveDigis + process.test + process.CutAs)
#process.schedule = cms.Schedule(process.p,process.endjob_step)
