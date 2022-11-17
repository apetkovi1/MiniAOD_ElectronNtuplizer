import CRABClient
from CRABClient.UserUtilities import config 

config = config()

#config.General.requestName = 'Bpark_MiniAOD_MC_reliso_UL' 
config.General.requestName = 'Bpark_MiniAOD_DATA_BPH4_runA_reliso_UL_v2' 
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py'


#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM'
#config.Data.inputDataset = '/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM'
config.Data.inputDataset = '/ParkingBPH4/Run2018A-UL2018_MiniAODv2-v1/MINIAOD' #'/ParkingBPH6/Run2018B-05May2019-v2/MINIAOD'                                      #'/ParkingBPH6/Run2018B-05May2019-v2/MINIAOD'
#config.Data.inputDataset = '/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIISummer19UL18MiniAODv2-TrkExtra_106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
#config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5 #
#config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
#config.Data.runRange = '275776-275782'
config.Data.publication = True
#config.Data.outputDatasetTag = 'Bpark_MiniAOD_MC_reliso_UL'
config.Data.outputDatasetTag = 'Bpark_MiniAOD_DATA_BPH4_runA_reliso_UL_v2'

config.Site.storageSite = 'T3_CH_CERNBOX'
