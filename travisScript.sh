#!/bin/bash

cd $1
#source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh

cmsrel CMSSW_9_3_3
cd CMSSW_9_3_3/src
cmsenv
cd -

g++ -g -O3 `root-config --glibs --cflags` fittingClass/mipFitsSiPM.C -o mipFitsSiPM
