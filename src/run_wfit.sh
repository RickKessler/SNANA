#!/usr/bin/env bash

clear
make ../bin/wfit.exe > tmp
#tail -n 10 tmp
../bin/wfit.exe /global/cscratch1/sd/kessler/PIPPIN_OUTPUT/PLASTICC_SPECTROZ_NEW/6_BIASCOR/BBC_SIMDATA/output/OUTPUT_BBCFIT/FITOPT001_MUOPT000.M0DIF  > log

tail -n 35 log
echo '--------------------------------------------------------------' 
wfit.exe /global/cscratch1/sd/kessler/PIPPIN_OUTPUT/PLASTICC_SPECTROZ_NEW/6_BIASCOR/BBC_SIMDATA/output/OUTPUT_BBCFIT/FITOPT001_MUOPT000.M0DIF > log2

tail -n 15 log2

