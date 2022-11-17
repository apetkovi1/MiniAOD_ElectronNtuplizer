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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", fileNames =
cms.untracked.vstring('/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/100000/041EF22D-69F5-914D-AD60-F2D1187B0842.root') #MC JPsi
#cms.untracked.vstring('/store/data/Run2018B/ParkingBPH6/MINIAOD/05May2019-v2/240000/005DC94B-063D-7245-8967-C9664533B57B.root') #data BParking
)

# /store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1_ext1-v1/40000/0076B640-F5D6-8545-B258-1ED605D2E791.root
# /store/data/Run2018B/ParkingBPH6/MINIAOD/05May2019-v2/240000/005DC94B-063D-7245-8967-C9664533B57B.root
# /store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/100000/041EF22D-69F5-914D-AD60-F2D1187B0842.root
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

