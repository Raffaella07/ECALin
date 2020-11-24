#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc700
cd /afs/cern.ch/work/r/ratramon/Bparking/CMSSW_10_2_15/src/
eval `scramv1 runtime -sh`
echo $PWD	

export ID=$(ls /eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct22_withB/ParkingBPH"$2"/crab_data_Run2018"$1"_part"$2"/)	
echo $ID
cd /afs/cern.ch/work/r/ratramon/ECAL_linearity/ECALinAnalysis/ECALlinearity/analysis

./Skim $1 $2 $3 $4 $ID
							
							
cd /afs/cern.ch/work/r/ratramon/ECAL_linearity/ECALinAnalysis/ECALlinearity/script
							
							
							
							
							
							
							
							
							
