from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'muMuGamma_ParkingPFCands_14Aug2023_1'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muMuGammaTree.py'

config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/ScoutingPFMonitor/Run2022F-v1/RAW'

# These values only make sense for processing data
#    Select input data based on a lumi mask
config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_US_MIT'
