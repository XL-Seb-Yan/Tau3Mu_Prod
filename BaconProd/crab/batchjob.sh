source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
cd /afs/cern.ch/user/x/xuyan/3MuonProj/CMSSW_8_0_27/src/
eval `scram runtime -sh`
cd BaconProd/crab/
cmsRun makingBaconPuppiMVAMets_MC.py
