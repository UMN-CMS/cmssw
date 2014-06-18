import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('tester',
                      bob = cms.double(10)
)
process.TFileService = cms.Service("TFileService",
                                          fileName = cms.string("/home/user1/lesko/work/CMSSW_7_1_0_pre5_HFTP/src/CalibCalorimetry/CaloTPG/hope.root")
                                          )
                                   
