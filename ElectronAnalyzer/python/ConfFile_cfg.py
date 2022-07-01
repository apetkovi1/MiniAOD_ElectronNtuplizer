import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process = cms.Process("Electrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v18')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", fileNames =
cms.untracked.vstring('file:00E6FFAF-1E0A-5049-9C12-FB3C7FFA466C.root') #MC JPsi
#cms.untracked.vstring('file:00496A25-08B6-FB4E-9681-D5FF4E1BE81F.root') #data BParking
)

setupEgammaPostRecoSeq(process,era='2018-Prompt')

process.electrons = cms.EDAnalyzer('ElectronAnalyzer',elecSrc = cms.untracked.InputTag("slimmedElectrons"),rhoSrc = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                                   pileupSrc = cms.untracked.InputTag("slimmedAddPileupInfo"))

#process.triggers = cms.EDAnalyzer("TriggerAnalyzer",
#    bits = cms.InputTag("TriggerResults","","HLT"),
#    prescales = cms.InputTag("patTrigger"),
#    objects = cms.InputTag("slimmedPatTrigger"),
#)

process.TFileService = cms.Service("TFileService", fileName=cms.string("TnP_tree.root"))

process.p = cms.Path(process.electrons*process.egammaPostRecoSeq)
