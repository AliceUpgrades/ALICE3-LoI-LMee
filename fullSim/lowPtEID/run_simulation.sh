#!/bin/bash
#
# This macro runs all steps in a strangeness tracking effort
#
# Step 1: create simulation via o2-sim
#
# Basic configurations
NEVENTS=50
NCPUS=1

DIR=`pwd`

# options to pass to every workflow
gloOpt=" -b --run --shm-segment-size 10000000000"

echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
o2-sim-serial -j ${NCPUS} --field 2 -e TGeant3 -n ${NEVENTS} -g pythia8hi -m TRK A3IP -o o2sim --trigger external --configKeyValues 'TriggerExternal.fileName=trigger_multiplicity.C;TriggerExternal.funcName=trigger_multiplicity(-1.5,1.5, 11)'
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
