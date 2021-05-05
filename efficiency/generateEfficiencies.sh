#! /usr/bin/env bash

runDelphes() {
  ### copy pythia8 configuration and adjust it
  if [[ $5 == "PbPb" ]]
  then
    cp ../pythia/pythia8.$5.cfg pythia8.$5.$1.cfg
    sleep 2
    echo "Main:numberOfEvents $2" >> pythia8.$5.$1.cfg
    echo "Random:seed = `expr 1001 \* $1`" >> pythia8.$5.$1.cfg
    echo "Random:setSeed on" >> pythia8.$5.$1.cfg
    echo "Beams:allowVertexSpread on " >> pythia8.$5.$1.cfg
    echo "Beams:sigmaTime 60." >> pythia8.$5.$1.cfg
  elif [[ $5 == "pp" ]]
  then
   cp ../pythia/pythia8.$5.default.cfg pythia8.$5.$1.cfg
  # cp ../pythia/pythia8.$5.cc.cfg pythia8.$5.cc.$1.cfg
   cp ../pythia/pythia8.$5.bb.cfg pythia8.$5.bb.$1.cfg
   sleep 2
   # write function to be more readable
   echo "Main:numberOfEvents $2" >> pythia8.$5.$1.cfg
   # echo "Main:numberOfEvents $3" >> pythia8.$5.cc.$1.cfg
   echo "Main:numberOfEvents $4" >> pythia8.$5.bb.$1.cfg
   # echo "Random:seed = $1$1$1" >> pythia8.$5.$1.cfg
   # echo "Random:seed = $1$1$1" >> pythia8.$5.cc.$1.cfg
   # echo "Random:seed = $1$1$1" >> pythia8.$5.bb.$1.cfg
   echo "Random:seed = `expr 1001 \* $1`" >> pythia8.$5.$1.cfg
   # echo "Random:seed = `expr 3002 \* $1`" >> pythia8.$5.cc.$1.cfg
   echo "Random:seed = `expr 5003 \* $1`" >> pythia8.$5.bb.$1.cfg
   echo "Beams:allowVertexSpread on " >> pythia8.$5.$1.cfg
   echo "Beams:allowVertexSpread on " >> pythia8.$5.bb.$1.cfg
   echo "Beams:sigmaTime 60." >> pythia8.$5.$1.cfg
   echo "Beams:sigmaTime 60." >> pythia8.$5.bb.$1.cfg
 else
   echo " !!! collisions System not available."
 fi



  # DelphesPythia8 propagate.tcl pythia8.$5.$1.cfg delphes.default.$1.root  &> delphes.default.$1.log &&
  # DelphesPythia8 propagate.tcl pythia8.$5.cc.$1.cfg delphes.cc.$1.root  &> delphes.cc.$1.log &&
  # DelphesPythia8 propagate.tcl pythia8.$5.bb.$1.cfg delphes.bb.$1.root  &> delphes.bb.$1.log &&
  DelphesPythia8 propagate.tcl pythia8.$5.$1.cfg delphes.$5.$1.root  &> delphes.$5.$1.log &&
  hadd -f delphes.$1.root delphes.*.$1.root && rm delphes.*.$1.root &&
  root -b -q -l "anaEEstudy.cxx(\"delphes.$1.root\", \"anaEEstudy.$1.root\")" &> anaEEstudy.$1.log
  # root -b -q -l "anaEEstudy.cxx(\"delphes.$1.root\", \"anaEEstudy.$1.root\")"
}
NJOBS=6        # number of max parallel runs
NRUNS=10        # number of runs

NEVENTS=10    # number of events in a run
NEVENTSCC=1000  # number of events in the charm sample
NEVENTSBB=1000  # number of events in the beauty sample

SYSTEM="PbPb"         # collisionSystem
# SYSTEM="pp"         # collisionSystem
# SCENARIO="default"     # detector setup
SCENARIO="Werner"     # detector setup
# BFIELD=2       # magnetic field  [kG]
BFIELD=5       # magnetic field  [kG]

SIGMAT=0.020      # time resolution [ns]
SIGMA0=0.200      # vertex time spread [ns]
BARRELRAD=100.    # barrel radius      [cm] (right now equal to TOF)
BARRELLEN=200.    # barrel half length [cm] (right now equal to TOF)
BARRELETA=1.443   # barrel max pseudorapidity
TAILLX=1.0        # tail on left    [q]
TAILRX=1.3        # tail on right   [q]
TOFRAD=100.       # TOF radius      [cm]
TOFLEN=200.       # TOF half length [cm]
RICHRAD=100.      # RICH radius      [cm]
RICHLEN=200.      # RICH half length [cm]

### calculate max eta from geometry
BARRELETA=`awk -v a=$TOFRAD -v b=$TOFLEN 'BEGIN {th=atan2(a,b)*0.5; sth=sin(th); cth=cos(th); print -log(sth/cth)}'`

#how many events are generated
echo " --- generating with scenario $SCENARIO setup"
if [[ $SYSTEM = "pp" ]]
  then
  ALLEVENTSINEL=$(expr $NEVENTS \* $NRUNS)
  ALLEVENTSCC=$(expr $NEVENTSCC \* $NRUNS)
  ALLEVENTSBB=$(expr $NEVENTSBB \* $NRUNS)
  echo " --- $ALLEVENTSINEL inelastic events"
  echo " --- $ALLEVENTSCC charm events"
  echo " --- $ALLEVENTSBB beauty events"
elif [[ $SYSTEM = "PbPb" ]]
 then
 ALLEVENTS=$(expr $NEVENTS \* $NRUNS)
 echo " --- $ALLEVENTS PbPb events"
fi
# copy stuff into running directory

echo " --- selected SYSTEM:   $SYSTEM"

# card
if [[ $BFIELD -eq 2 ]]
then
  cp ../delphes/cards/propagate.2kG.tails.tcl propagate.tcl
  echo " --- selected B-Field:  0.${BFIELD}T"
elif [[ $BFIELD -eq 5 ]]
then
  cp ../delphes/cards/propagate.5kG.tails.tcl propagate.tcl
  echo " --- selected B-Field:  0.${BFIELD}T"
else echo " !!! BFIELD not available, check vstd::coutaribale"
fi

echo " --- selected SCENARIO: $SCENARIO"

# code
cp ./macros/anaEEstudy.cxx anaEEstudy.cxx

# LUTS
if [[ $SCENARIO = "Werner" ]] && [[ $BFIELD -eq 2 ]]
then
  cp ../LUTs/lutCovm.werner.rmin100.2kG/lutCovm.el.werner.rmin100.2kG.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.werner.rmin100.2kG/lutCovm.mu.werner.rmin100.2kG.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.werner.rmin100.2kG/lutCovm.pi.werner.rmin100.2kG.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.werner.rmin100.2kG/lutCovm.ka.werner.rmin100.2kG.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.werner.rmin100.2kG/lutCovm.pr.werner.rmin100.2kG.dat lutCovm.pr.dat
elif [[ $SCENARIO = "Werner" ]] && [[ $BFIELD -eq 5 ]]
then
  cp ../LUTs/lutCovm.werner.rmin100.5kG/lutCovm.el.werner.rmin100.5kG.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.werner.rmin100.5kG/lutCovm.mu.werner.rmin100.5kG.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.werner.rmin100.5kG/lutCovm.pi.werner.rmin100.5kG.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.werner.rmin100.5kG/lutCovm.ka.werner.rmin100.5kG.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.werner.rmin100.5kG/lutCovm.pr.werner.rmin100.5kG.dat lutCovm.pr.dat
elif [[ $SCENARIO = "default" ]] && [[ $BFIELD -eq 2 ]]
then
  cp ../LUTs/lutCovm.2kG.100cm.default/lutCovm.el.2kG.100cm.default.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.2kG.100cm.default/lutCovm.mu.2kG.100cm.default.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.2kG.100cm.default/lutCovm.pi.2kG.100cm.default.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.2kG.100cm.default/lutCovm.ka.2kG.100cm.default.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.2kG.100cm.default/lu      # time resolution [ns]
      # vertex time spread [ns]
    # barrel radius      [cm] (right now equal to TOF)
    # barrel half length [cm] (right now equal to TOF)
   # TOF max pseudorapidity
        # tail on left    [q]
        # tail on right   [q]
       # TOF radius      [cm]
       # TOF half length [cm]
      # RICH radius      [cm]
      # RICH half length [cm]tCovm.pr.2kG.100cm.default.dat lutCovm.pr.dat
elif [[ $SCENARIO = "default" ]] && [[ $BFIELD -eq 5 ]]
then
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.el.5kG.100cm.default.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.mu.5kG.100cm.default.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.pi.5kG.100cm.default.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.ka.5kG.100cm.default.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.pr.5kG.100cm.default.dat lutCovm.pr.dat
else
  echo "!!! check SCENARIO and BFIELD variables"
fi

# make sure B field is set right
sed -i -e "s/set Bz .*$/set Bz ${BFIELD}e\-1/" propagate.tcl
sed -i -e "s/double Bz .*$/double Bz = ${BFIELD}e\-1;/" anaEEstudy.cxx
### set TOF radius
sed -i -e "s/set barrel_Radius .*$/set barrel_Radius ${BARRELRAD}e\-2/" propagate.tcl
sed -i -e "s/double tof_radius = .*$/double tof_radius = ${TOFRAD}\;/" anaEEstudy.cxx
### set TOF length
sed -i -e "s/set barrel_HalfLength .*$/set barrel_HalfLength ${BARRELLEN}e\-2/" propagate.tcl
sed -i -e "s/double tof_length = .*$/double tof_length = ${TOFLEN}\;/" anaEEstudy.cxx
### set TOF acceptance
sed -i -e "s/set barrel_Acceptance .*$/set barrel_Acceptance \{ 0.0 + 1.0 * fabs(eta) < ${BARRELETA} \}/" propagate.tcl
### set TOF time resolution and tails
sed -i -e "s/set barrel_TimeResolution .*$/set barrel_TimeResolution ${SIGMAT}e\-9/" propagate.tcl
sed -i -e "s/set barrel_TailRight .*$/set barrel_TailRight ${TAILRX}/" propagate.tcl
sed -i -e "s/set barrel_TailLeft  .*$/set barrel_TailLeft ${TAILLX}/" propagate.tcl
sed -i -e "s/double tof_sigmat = .*$/double tof_sigmat = ${SIGMAT}\;/" anaEEstudy.cxx
sed -i -e "s/double tof_sigma0 = .*$/double tof_sigma0 = ${SIGMA0}\;/" anaEEstudy.cxx
### set RICH radius
sed -i -e "s/double rich_radius = .*$/double rich_radius = ${RICHRAD}\;/" anaEEstudy.cxx
### set RICH length
sed -i -e "s/double rich_length = .*$/double rich_length = ${RICHLEN}\;/" anaEEstudy.cxx


# adapt pt cuts corresponding to selected B-field
if [[ $BFIELD -eq 2 ]]
then
  PTACC=4 # pt = 40 MeV/c acceptance cut for low B field
elif [[ $BFIELD -eq 5 ]]
then
  PTACC=8 # pt = 40 MeV/c acceptance cut for low B field
fi
sed -i -e "s/double PtCut .*$/double PtCut = ${PTACC}e\-2;/" anaEEstudy.cxx
echo " --- using Pt Acceptance cut:  0.0${PTACC}GeV/c"

#load the right LUTs in the anaEEstudy.cxx

cp ../preshower/macros/preparePreshowerEff.C .
root -l -b -q preparePreshowerEff.C
echo "Finished   preparePreshowerEff.C"
echo ""


echo " ----------------------------------"
echo " List of set detector parameters: "
echo " SIGMAT      = ${SIGMAT}      # time resolution [ns]"
echo " SIGMA0      = ${SIGMA0}      # vertex time spread [ns]"
echo " BARRELRAD   = ${BARRELRAD}       # barrel radius      [cm] (right now equal to TOF)"
echo " BARRELLEN   = ${BARRELLEN}       # barrel half length [cm] (right now equal to TOF)"
echo " BARRELETA   = ${BARRELETA}    # barrel max pseudorapidity"
echo " TAILLX      = ${TAILLX}        # tail on left    [q]"
echo " TAILRX      = ${TAILRX}        # tail on right   [q]"
echo " TOFRAD      = ${TOFRAD}       # TOF radius      [cm]"
echo " TOFLEN      = ${TOFLEN}       # TOF half length [cm]"
echo " RICHRAD     = ${RICHRAD}       # RICH radius      [cm]"
echo " RICHLEN     = ${RICHLEN}       # RICH half length [cm]"
echo " ----------------------------------"
echo ""

### loop over runs
echo " Run Delphes... "
for I in $(seq 1 $NRUNS); do

    ### wait for a free slot
    while [ $(ls .running.* 2>/dev/null | wc -l) -ge $NJOBS ]; do
      echo " --- waiting for a free slot"
      sleep 10
    done

    ### book the slot
    echo " --- starting run $I"
    touch .running.$I

    runDelphes $I $NEVENTS $NEVENTSCC $NEVENTSBB $SYSTEM &&
    (rm -rf delphes.$I.root && rm -rf .running.$I && echo " --- run $I completed") ||
    (rm -rf delphes.$I.root && rm -rf .running.$I && echo " --- run $I crashed") &

done
echo " Finished running Delphes "


### merge runs when all done
wait
hadd -f anaEEstudy.${SYSTEM}.${SCENARIO}.B=0.${BFIELD}_$(expr $NEVENTS \* $NRUNS)events.root anaEEstudy.*.root && #rm -rf anaEEstudy.*.root &&
### run analysis scrip on new file
cp ./macros/anaPlots.cxx anaPlots.cxx &&
root -l -b -q "anaPlots.cxx(\"anaEEstudy.${SYSTEM}.${SCENARIO}.B=0.${BFIELD}_$(expr $NEVENTS \* $NRUNS)events.root\")" &&
rm anaPlots.cxx
mv anaEEstudy.${SYSTEM}.${SCENARIO}.B=0.${BFIELD}_$(expr $NEVENTS \* $NRUNS)events.root ./data/prod
### clean up
rm lutCovm*
rm propagate.tcl
rm *.root
rm *.log
rm *.cfg
rm anaEEstudy.cxx
rm preparePreshowerEff.C
# ./cleanup.sh
