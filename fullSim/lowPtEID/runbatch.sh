#!/bin/bash
echo "CPU #${1} starting"
mkdir job_${1}
cp run_simulation.sh job_${1}
cp RecoDecay.h job_${1}
cp trigger_multiplicity.C job_${1}
cp run_strangenesstracking_01.C job_${1}
cd job_${1}
source run_simulation.sh
root.exe -q -b run_strangenesstracking_01.C+
if [ ! -f histooutput.root ]; then
        echo "Failed, trying again"
        root.exe -q -b run_strangenesstracking_01.C+
fi
if [ ! -f histooutput.root ]; then
        echo "Failed yet again. Third time's the charm"
        root.exe -q -b run_strangenesstracking_01.C+
fi
if [ ! -f histooutput.root ]; then
        echo "Failed yet again. Third time's the charm"
        root.exe -q -b run_strangenesstracking_01.C+
fi
if [ ! -f histooutput.root ]; then
        echo "Forget it. Clean tree output up, this is bogus, sorry"
        rm treeoutput.root
fi
echo "Doing master cleanup"
rm o2sim_Kine.root
rm o2sim_HitsTRK.root
rm o2sim.root
rm run_strangenesstracking_01_C.d
rm run_strangenesstracking_01_C.so
rm o2sim_par.root
rm o2sim_geometry.root
rm o2sim_MCHeader.root
rm o2sim_configuration.ini
rm run_strangenesstracking_01_C_ACLiC_dict_rdict.pcm
rm o2sim_grp.root
rm MCStepLoggerVolMap.dat
rm o2sim_proc-cut.dat
rm MCStepLoggerSenVol.dat
rm gphysi.dat
echo "Done! Enjoy!"
