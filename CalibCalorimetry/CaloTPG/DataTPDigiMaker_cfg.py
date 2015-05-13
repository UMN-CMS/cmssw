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
    TempFN=['file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/02FB3077-53A6-E111-934B-BCAEC5329700.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/0456F6A3-2CA6-E111-A084-5404A63886B6.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/06C977BC-5EA6-E111-B4B7-0025901D5C86.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/06EBF815-2BA6-E111-B535-BCAEC518FF54.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/0EAA7493-1FA6-E111-A57F-5404A63886A2.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/140F5365-1DA6-E111-9AC9-003048F024FE.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/14DCCFE3-59A6-E111-8576-0025901D5C86.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/18316693-1FA6-E111-9C7B-BCAEC518FF50.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/28ECC3D6-30A6-E111-82B3-003048F117EC.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/2A34EBB3-3CA6-E111-8E6F-0025901D626C.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/3039EE8E-3FA6-E111-8C4C-5404A63886B4.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/344DC8ED-1BA6-E111-BC0C-5404A63886A0.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/42303B27-65A6-E111-950F-5404A63886BE.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/447017D0-3EA6-E111-A972-003048D2BE08.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/44EB6113-18A6-E111-8205-5404A63886AE.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/4CF54C9F-5CA6-E111-A59F-BCAEC5364C42.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/5275E88D-1AA6-E111-9E2A-001D09F29524.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/52A3E7CF-23A6-E111-8344-0025B32036D2.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/686E41F7-34A6-E111-95A6-003048F1C420.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/6A87F431-4AA6-E111-A93E-E0CB4E4408E3.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/76167738-71A6-E111-B7C7-BCAEC532971C.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/78171AA2-25A6-E111-A8F5-001D09F23F2A.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/7819B206-63A6-E111-81AD-003048678098.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/7C2DE175-16A6-E111-BF07-003048D2C174.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/7C318BE8-4FA6-E111-8071-BCAEC5329719.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/7CD74957-3BA6-E111-A2A8-5404A63886B0.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/8224362D-30A6-E111-95B3-002481E94C7E.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/881A1431-2FA6-E111-9036-BCAEC5364C42.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/8A3FD822-43A6-E111-9245-003048D2C0F2.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/92A95962-4CA6-E111-BD4B-BCAEC518FF8F.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/9AE494B7-14A6-E111-9A2A-003048F024DA.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/A609B9AD-2DA6-E111-AD0F-5404A63886A2.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/A89F0B90-15A6-E111-862D-001D09F23A20.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/B0B9BDAB-57A6-E111-AEB6-0025901D6268.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/B4374C02-19A6-E111-B8C2-002481E0D83E.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/B8A28509-37A6-E111-88EF-5404A638869E.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/B8CAB517-12A6-E111-9934-E0CB4E55365D.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/BC11AA2E-39A6-E111-B74C-BCAEC5364CED.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/BE3CB835-56A6-E111-8672-BCAEC518FF76.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/CA629788-33A6-E111-9B03-003048F024DE.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/CE83E9B0-20A6-E111-90F2-BCAEC518FF30.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/D4702F14-48A6-E111-A1CA-00215AEDFD74.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/D848E8B2-46A6-E111-A489-0025901D623C.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/D8A9AD6E-32A6-E111-9BA8-0025901D5D7E.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/D8BAD72A-27A6-E111-9F6C-003048CF99BA.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/DE875A32-1AA6-E111-9EDB-0025B324400C.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/DEDD4FAD-35A6-E111-A8F6-0025901D5DF4.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/E0AA8E53-22A6-E111-9ED2-00215AEDFD98.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/E2603AAF-44A6-E111-8241-003048F117B6.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/E279E4B7-28A6-E111-BB00-5404A638869C.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/E2B9C710-4DA6-E111-80F8-BCAEC518FF80.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/E464B525-52A6-E111-A94C-0025901D5D80.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/EC2D6375-38A6-E111-BB91-0025901D5E10.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/F2E213FA-60A6-E111-9498-001D09F27003.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/F4BB01C1-12A6-E111-A8F4-001D09F25267.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/FE1AE202-41A6-E111-935C-5404A63886C1.root',
    'file:/hdfs/cms/phedex/store/data/Run2012B/MinimumBias/RAW/v1/000/194/912/FEBD36C8-6FA6-E111-B0D6-001D09F295A1.root']

              
      
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



