#!/bin/bash
#eval $(alienv printenv VO_ALICE@DelphesO2::v20210817-1)
TMPDIR=$PWD
PBS_O_WORKDIR=$PWD
cd ${TMPDIR}
echo ${TMPDIR}
echo ${PBS_O_WORKDIR}

NEVENTS=5
NCPUS=1

# options to pass to every workflow
gloOpt=" -b --run --shm-segment-size 10000000000"
export VMCWORKDIR=${O2_ROOT}/share

outdir=${PBS_O_WORKDIR}/out

mkdir -p ${outdir}


#cp ${PBS_O_WORKDIR}/run_stratrack.sh ${TMPDIR}
#cp ${PBS_O_WORKDIR}/config_custom_Pythia.cfg ${TMPDIR}
#cp ${PBS_O_WORKDIR}/configPythia.ini ${TMPDIR}
#cp ${PBS_O_WORKDIR}/generator_pythia8_gun.C ${TMPDIR}
#cp ${PBS_O_WORKDIR}/RecoDecay.h ${TMPDIR}
#cp ${PBS_O_WORKDIR}/pythia8_hi.cmnd ${TMPDIR}
#cp ${PBS_O_WORKDIR}/pp13.cmnd ${TMPDIR}
#cp ${PBS_O_WORKDIR}/run_strangenesstracking_xicc.C ${TMPDIR}
#cp ${PBS_O_WORKDIR}/trigger_multiplicity.C ${TMPDIR}

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

rm run_strangenesstracking_xicc.C~
rm run_strangenesstracking_xicc_C_ACLiC_dict_rdict.pcm

mv vertexer.qa.root ${outdir}
mv xicc.treeoutput.root ${outdir}
mv xicc.histooutput.root ${outdir}
echo "Anything left?" 
ls
echo "Done! Enjoy!"
