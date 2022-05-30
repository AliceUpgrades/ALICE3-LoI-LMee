#!/bin/bash
#
# This macro runs all steps in a strangeness tracking effort
#
# Step 1: create simulation via o2-sim
#
# Basic configurations
NEVENTS=1000
NCPUS=1

DIR=`pwd`

# options to pass to every workflow
gloOpt=" -b --run --shm-segment-size 10000000000"

echo "CPU #${1} starting"
mkdir job_${1}
cp run_stratrack.sh job_${1}
cp config_custom_Pythia.cfg job_${1}
cp configPythia.ini job_${1}
cp generator_pythia8_PbPb.C job_${1}
cp RecoDecay.h job_${1}
cp pythia8_hi.cmnd job_${1}
cp run_strangenesstracking_xi.C job_${1}
cp findable.C job_${1}
cd job_${1}

#Roberto's functions
## create the transport.C macro
cat <<EOF > transport.C
o2::data::Stack::TransportFcn
transport()
{
  return [](const TParticle& p, const std::vector<TParticle>& particles) -> bool {
           auto eta = p.Eta();
           if (std::fabs(eta) > 1.5) return false;
           auto pdg = std::abs(p.GetPdgCode());
           if (pdg == 310) return true;
           if (pdg == 3122) return true;
           if (pdg == 3312) return true;
           if (pdg == 3334) return true;
           return false;
   };
}
EOF

##run the simulation
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
time o2-sim-serial -j ${NCPUS} --field 20U -e TGeant3 -n ${NEVENTS} -g external --configFile configPythia.ini -m A3IP TRK --configKeyValues "Diamond.width[2]=6.;Diamond.width[0]=0.005;Diamond.width[1]=0.005;" -o o2sim
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
#root.exe -q -b run_strangenesstracking_xi.C+
#if [ ! -f xi.histooutput.root ]; then
#        echo "Failed, trying again"
#        root.exe -q -b run_strangenesstracking_xi.C+
#fi
#if [ ! -f xi.histooutput.root ]; then
#        echo "Failed yet again. Third time's the charm"
#        root.exe -q -b run_strangenesstracking_xi.C+
#fi
#if [ ! -f xi.histooutput.root ]; then
#        echo "Failed yet again. Third time's the charm"
#        root.exe -q -b run_strangenesstracking_xi.C+
#fi
#if [ ! -f xi.histooutput.root ]; then
#        echo "Forget it. Clean tree output up, this is bogus, sorry"
#        rm xi.treeoutput.root
#fi

root.exe -q -b findable.C+
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
echo "Doing master cleanup"
#rm o2sim_Kine.root
#rm o2sim_HitsTRK.root
rm o2sim.root
rm run_strangenesstracking_xi_C.d
rm run_strangenesstracking_xi_C.so
rm o2sim_par.root
#rm o2sim_geometry.root
rm o2sim_MCHeader.root
rm o2sim_configuration.ini
rm run_strangenesstracking_xi_C_ACLiC_dict_rdict.pcm
#rm o2sim_grp.root
rm MCStepLoggerVolMap.dat
rm o2sim_proc-cut.dat
rm MCStepLoggerSenVol.dat
rm gphysi.dat
echo "Done! Enjoy!"
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
