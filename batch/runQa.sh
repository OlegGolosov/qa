#!/bin/bash

date
source ${rootConfig}

jobId=$SLURM_ARRAY_TASK_ID
dataFile=$(sed "${jobId}q;d" $fileList)

qaFile=$qaDir/$(basename $dataFile)
qaFile=${qaFile/.tree.root/.qa.root}
qaFile=${qaFile/.list/.qa.root}

echo Processing file: $dataFile 
echo QA file: $qaFile 

cd $qaDir/src
. env.sh
time root -b -l -q qa.C"(\"$dataFile\",\"$qaFile\",\"$eventSelection\")"

exit $?
