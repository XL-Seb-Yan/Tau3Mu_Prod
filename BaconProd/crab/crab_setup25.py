from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'data7'
config.General.workArea =  'runf2'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'makingBaconPuppiMVAMets_MC.py'
#config.JobType.psetName = 'makingGenBacon.py'
config.JobType.psetName = 'makingBaconPuppiMVAMets_DATA.py'
config.JobType.outputFiles = ['Output.root']
config.section_("Data")
#config.Data.inputDataset = '/SingleMuon/Run2016B-PromptReco-v1/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2016B-PromptReco-v2/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2016C-PromptReco-v2/AOD'
#config.Data.inputDataset = '/DsToTauNeutrinoToMuMuMuNeutrino_TuneCUEP8M1_13TeV-pythia8-evtGen/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'
#config.Data.inputDataset = '/DsToTauToMuMuMu_MuFilter_TuneCUEP8M1_13TeV-pythia8/RunIISummer16DR80-FlatPU28to62HcalNZSRAWAODSIM_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/RAWAODSIM'
#config.Data.inputDataset = '/DoubleMuon/Run2016H-23Sep2016-v1/AOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016H-PromptReco-v3/AOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016G-18Apr2017-v1/AOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016B-18Apr2017-v1/AOD'
#config.Data.inputDataset = '/SingleElectron/Run2016B-PromptReco-v2/AOD'
#config.Data.inputDataset = '/SingleElectron/Run2016B-PromptReco-v1/AOD'
#config.Data.inputDataset = '/DoubleEG/Run2016B-PromptReco-v2/AOD'
#config.Data.inputDataset = '/SingleElectron/Run2016B-PromptReco-v1/AOD'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
#config.Data.inputDataset = '/VBF_HToInvisible_M125_13TeV_powheg_pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'
#config.Data.inputDataset = '/SingleElectron/Run2016H-PromptReco-v1/AOD'
config.Data.inputDataset = '/SingleElectron/Run2016F-23Sep2016-v1/AOD'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/AODSIM'

#config.JobType.inputFiles  = ['Summer15_25nsV6_DATA.db']
#config.JobType.inputFiles  = ['Summer15_25nsV6_MC.db']
#config.Data.ignoreLocality = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.lumiMask = 'rereco2.json'
config.Data.unitsPerJob = 7
#config.Data.outLFNDirBase = '/store/group/comm_luminosity/SingleLepton'
config.Data.outLFNDirBase = '/store/user/arapyan/Run2' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'Samples'
#config.Data.publishDataName = 'Samples'

config.section_("Site")
config.Site.storageSite = 'T3_US_MIT'
