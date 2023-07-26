#!/bin/bash

dataDir=$1
listDir=$2
nLines=$3

mkdir -p $listDir

$dataDir/*root > $listDir/all.list
split -l $nLines --numeric-suffixes=1 --additional-suffix=.list $listDir/all.list $listDir/
rm $listDir/all.list
