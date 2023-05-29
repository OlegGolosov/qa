#!/bin/bash

system=pbpb
pbeam=13
production=dcmqgsm_smm
export eventSelection=goodMcEvent
postfix=
dataDir=/mnt/pool/nica/7/ovgol/na61_data/${system}/${pbeam}agev/${production}
export qaDir=/mnt/pool/nica/7/ovgol/na61_qa/${system}/${pbeam}agev/${production}/${eventSelection}/${postfix}
logDir=/mnt/pool/nica/7/ovgol/log
dataPattern=all.list
#dataPattern=*.tree.root

export rootConfig=/mnt/pool/nica/7/parfenovpeter/Soft/root-6-26/install/bin/thisroot.sh
qaSrc=/home/ovgol/soft/qa
executable=./runQa.sh

mkdir -pv $logDir

mkdir -pv $qaDir/src
cp -v $0 $qaDir/src
cp -v $executable $qaDir/src
executable=$qaDir/src/$(basename $executable)
qaSources=(\
  env.sh\
  qa.C\
  qaUtils.h\
  histsKIT.h
  histsDT.h
  makeDF.h\
  vtxPurity_pbpb13.root\
  pid_pbpb13_graphical.root\
  pbpb13.pid.root\
  centrality_pbpb13_hE.root\
  cut_dEdx_m2_eneg.C\       
  cut_dEdx_m2_epos.C\     
  cut_dEdx_m2_kaonneg.C\  
  cut_dEdx_m2_kaonpos.C\
  cut_dEdx_m2_pionneg.C\
  cut_dEdx_m2_pionpos.C\
  cut_dEdx_m2_proton.C\
  cut_dEdx_m2_deuteron.C\
)
for f in ${qaSources[*]};do
  [ -e $f ] || cp -v $qaSrc/$f $qaDir/src
done

export fileList=${qaDir}/src/list
ls $dataDir/$dataPattern > $fileList
nFiles=$(cat $fileList|wc -w)

#sbatch -p cpu -a 1-$nFiles -o ${logDir}/qa_%A_%a.out $executable

export SLURM_ARRAY_TASK_ID=1
$executable
