import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/data/whybee0b/user/aevans/DataDIGIoutputRun200091.root'
    )
)

process.myProducerLabel = cms.EDProducer('PrimCopier'
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root'),
  SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring( 'keep *_*_*_OWNPARTICLES')
    
)


  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.p,process.e)
