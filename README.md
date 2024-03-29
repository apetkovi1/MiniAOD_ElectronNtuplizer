# MiniAOD_ElectronNtuplizer
## Description
This page contains instructions and examples of a way to extract information for Electron Tag and Probe from MiniAOD files

## Usage instructions

* Set up your area:
```
cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src
cmsenv
```
* Get the Egamma tool for extraction of MVA V2 value:
```
git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data -rf
git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
```
* Get the MiniAOD_ElectronNtuplizer:
```
git clone https://github.com/apetkovi1/MiniAOD_ElectronNtuplizer.git
```
* Compile everything and run (change the input file in ConfFile_cfg.py), the resulting ntuple should be stored in Tnp_tree.root
```
scram b
cd ElectronAnalyzer/python
cmsRun ConfFile_cfg.py
```

* If the Configuration files works fine and you get an output -> You can send CRAB jobs to get the result from the whole DATASET (all files in dataset). After running the submit, you'll be able to check on the progress of the jobs with line provided in terminal output.
```
voms-proxy-init --voms cms --valid 96:0
crab submit -c CrabConfig_mini.py
```
