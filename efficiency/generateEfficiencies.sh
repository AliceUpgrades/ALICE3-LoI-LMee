#! /usr/bin/env bash

runDelphes() {
  ### copy pythia8 configuration and adjust it
  if [[ $5 == "PbPb" ]]
  then
    cp ../pythia/pythia8.$5.cfg pythia8.$5.$1.cfg
    sleep 2
    echo "Main:numberOfEvents $2" >> pythia8.$5.$1.cfg
    echo "Random:seed = `expr 1001 \* $1`" >> pythia8.$5.$1.cfg
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
NJOBS=7        # number of max parallel runs
NRUNS=1        # number of runs

NEVENTS=100    # number of events in a run
NEVENTSCC=1000  # number of events in the charm sample
NEVENTSBB=1000  # number of events in the beauty sample

# SYSTEM="PbPb"         # collisionSystem
SYSTEM="pp"         # collisionSystem
# SCENARIO="default"     # detector setup
SCENARIO="Werner"     # detector setup
# BFIELD=2       # magnetic field  [kG]
BFIELD=5       # magnetic field  [kG]

#how many events are generated
ALLEVENTSINEL=$(expr $NEVENTS \* $NRUNS)
ALLEVENTSCC=$(expr $NEVENTSCC \* $NRUNS)
ALLEVENTSBB=$(expr $NEVENTSBB \* $NRUNS)
echo " --- generating with scenario $SCENARIO setup"
echo " --- $ALLEVENTSINEL inelastic events"
echo " --- $ALLEVENTSCC charm events"
echo " --- $ALLEVENTSBB beauty events"
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
  cp ../LUTs/lutCovm.2kG.100cm.default/lutCovm.pr.2kG.100cm.default.dat lutCovm.pr.dat
elif [[ $SCENARIO = "default" ]] && [[ $BFIELD -eq 5 ]]
then
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.el.5kG.100cm.default.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.mu.5kG.100cm.default.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.pi.5kG.100cm.default.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.ka.5kG.100cm.default.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.5kG.100cm.default/lutCovm.pr.5kG.100cm.default.dat lutCovm.pr.dat
else cp ./macros/anaPlots.cxx anaPlots.cxx &&
 root -l -b -q "anaPlots.cxx(\"anaEEstudy.${SYSTEM}.${SCENARIO}.B=0.${BFIELD}_$(expr $NEVENTS \* $NRUNS)events.root\", $SCENARIO)" &&

  echo "!!! check SCENARIO and BFIELD variables"
fi

# make sure B field is set right
sed -i -e "s/set Bz .*$/set Bz ${BFIELD}e\-1/" propagate.tcl
sed -i -e "s/double Bz .*$/double Bz = ${BFIELD}e\-1;/" anaEEstudy.cxx


#load the right LUTs in the anaEEstudy.cxx

### loop over runs
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
# ./cleanup.sh
