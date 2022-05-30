#!/bin/bash
eval $(alienv printenv VO_ALICE@DelphesO2::v20210817-1)

cd ${TMPDIR}
echo ${TMPDIR}
echo ${PBS_O_WORKDIR}

NEVENTS=500
NCPUS=1

# options to pass to every workflow
gloOpt=" -b --run --shm-segment-size 10000000000"
export VMCWORKDIR=${O2_ROOT}/share

outdir=${PBS_O_WORKDIR}/job_${PBS_ARRAYID}

mkdir -p ${outdir}

cp ${PBS_O_WORKDIR}/run_stratrack.sh ${TMPDIR}
cp ${PBS_O_WORKDIR}/config_custom_Pythia.cfg ${TMPDIR}
cp ${PBS_O_WORKDIR}/configPythia.ini ${TMPDIR}
cp ${PBS_O_WORKDIR}/generator_pythia8_gun.C ${TMPDIR}
cp ${PBS_O_WORKDIR}/RecoDecay.h ${TMPDIR}
cp ${PBS_O_WORKDIR}/pp13.cmnd ${TMPDIR}
cp ${PBS_O_WORKDIR}/run_strangenesstracking_xicc.C ${TMPDIR}

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

echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
time o2-sim-serial -j ${NCPUS} --field 20U -e TGeant3 -n ${NEVENTS} -g external --configFile configPythia.ini -m A3IP TRK -o o2sim --configKeyValues "Diamond.width[2]=6.;Diamond.width[0]=0.005;Diamond.width[1]=0.005;"
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
root.exe -q -b "run_strangenesstracking_xicc.C+(\"./\", \"${PBS_O_WORKDIR}/\")"
if [ ! -f xicc.histooutput.root ]; then
        echo "Failed, trying again"
	root.exe -q -b "run_strangenesstracking_xicc.C+(\"./\", \"${PBS_O_WORKDIR}/\")"
fi
if [ ! -f xicc.histooutput.root ]; then
        echo "Failed yet again. Third time's the charm"
	root.exe -q -b "run_strangenesstracking_xicc.C+(\"./\", \"${PBS_O_WORKDIR}/\")"
fi
if [ ! -f xicc.histooutput.root ]; then
        echo "Failed yet again. Third time's the charm"
	root.exe -q -b "run_strangenesstracking_xicc.C+(\"./\", \"${PBS_O_WORKDIR}/\")"
fi
if [ ! -f xicc.histooutput.root ]; then
        echo "Forget it. Clean tree output up, this is bogus, sorry"
        rm xicc.treeoutput.root
fi
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"


echo "Doing master cleanup"

rm MCStepLoggerSenVol.dat
rm MCStepLoggerVolMap.dat
rm RecoDecay.h
rm configPythia.ini
rm config_custom_Pythia.cfg
rm generator_pythia8_gun.C
rm gphysi.dat
rm o2sim.root
rm o2sim_HitsTRK.root
rm o2sim_Kine.root
rm o2sim_MCHeader.root
rm o2sim_configuration.ini
rm o2sim_geometry.root
rm o2sim_grp.root
rm o2sim_par.root
rm o2sim_proc-cut.dat
rm pp13.cmnd
rm run_strangenesstracking_xicc.C
rm run_strangenesstracking_xicc_C.d
rm run_strangenesstracking_xicc_C.so
rm run_strangenesstracking_xicc_C_ACLiC_dict_rdict.pcm
rm transport.C

mv xicc.treeoutput.root ${outdir}
mv vertexer.qa.root ${outdir}
mv xicc.histooutput.root ${outdir}

echo "Anything left?" 
ls
echo "Done! Enjoy!"
