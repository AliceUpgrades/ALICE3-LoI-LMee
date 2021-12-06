#! /usr/bin/env bash

runDelphes() {
  ### copy pythia8 configuration and adjust it
  if [[ $5 == "PbPb" ]]
  then
    cp ../pythia/pythia8_$5.cfg pythia8_$5.$1.cfg
    # cp ../pythia/pythia8_${5}_502TeV.cfg pythia8_$5.$1.cfg
    sleep 2
    echo "Main:numberOfEvents $2" >> pythia8_$5.$1.cfg
    echo "Random:setSeed on" >> pythia8_$5.$1.cfg
    echo "Random:seed = `expr 1001 \* $1`" >> pythia8_$5.$1.cfg
    # echo "Beams:allowVertexSpread on " >> pythia8_$5.$1.cfg
    # echo "Beams:sigmaTime 60." >> pythia8_$5.$1.cfg
  elif [[ $5 == "pp" ]]
  then
   cp ../pythia/pythia8_${5}_default.cfg pythia8_$5.$1.cfg
  # cp ../pythia/pythia8_${5}_cc.cfg pythia8_${5}_cc.$1.cfg
   cp ../pythia/pythia8_${5}_bb.cfg pythia8_${5}_bb.$1.cfg
   sleep 2
   # write function to be more readable
   echo "Main:numberOfEvents $2" >> pythia8_$5.$1.cfg
   # echo "Main:numberOfEvents $3" >> pythia8_${5}_cc.$1.cfg
   echo "Main:numberOfEvents $4" >> pythia8_${5}_bb.$1.cfg
   # echo "Random:seed = $1$1$1" >> pythia8_$5.$1.cfg
   # echo "Random:seed = $1$1$1" >> pythia8_$5.cc.$1.cfg
   # echo "Random:seed = $1$1$1" >> pythia8_$5.bb.$1.cfg
   echo "Random:setSeed on" >> pythia8_$5.$1.cfg
   # echo "Random:setSeed on" >> pythia8_${5}_cc.$1.cfg
   echo "Random:setSeed on" >> pythia8_${5}_bb.$1.cfg
   echo "Random:seed = `expr 1001 \* $1`" >> pythia8_$5.$1.cfg
   # echo "Random:seed = `expr 3002 \* $1`" >> pythia8_${5}_cc.$1.cfg
   echo "Random:seed = `expr 5003 \* $1`" >> pythia8_${5}_bb.$1.cfg
   # echo "Beams:allowVertexSpread on " >> pythia8_$5.$1.cfg
   # echo "Beams:allowVertexSpread on " >> pythia8_$5.bb.$1.cfg
   # echo "Beams:sigmaTime 60." >> pythia8_$5.$1.cfg
   # echo "Beams:sigmaTime 60." >> pythia8_$5.bb.$1.cfg
 else
   echo " !!! collisions System not available."
 fi



  # DelphesPythia8 propagate.tcl pythia8_$5.$1.cfg delphes.default.$1.root  &> delphes.default.$1.log &&
  # DelphesPythia8 propagate.tcl pythia8_$5.cc.$1.cfg delphes.cc.$1.root  &> delphes.cc.$1.log &&
  # DelphesPythia8 propagate.tcl pythia8_$5.bb.$1.cfg delphes.bb.$1.root  &> delphes.bb.$1.log &&
  DelphesPythia8 propagate.tcl pythia8_$5.$1.cfg delphes.$5.$1.root  &> delphes.$5.$1.log &&
  hadd -f delphes.$1.root delphes.*.$1.root && rm delphes.*.$1.root &&
  root -b -q -l "anaEEstudy.cxx(\"delphes.$1.root\", \"anaEEstudy.$1.root\")" &> anaEEstudy.$1.log
    # root -b -q -l "anaEEstudy.cxx(\"delphes.$1.root\", \"anaEEstudy.$1.root\")"
}


NJOBS=100        # number of max parallel runs                   50
NRUNS=100        # number of runs                               100

NEVENTS=1000    # number of events in a run                     1000
NEVENTSCC=1000  # number of events in the charm sample
NEVENTSBB=1000  # number of events in the beauty sample

SYSTEM="PbPb"            # collisionSystem
# SYSTEM="pp"            # collisionSystem
# SCENARIO="default"     # detector setup
SEMICENTRAL=true         # semicentral = false => 0-10%,  semicentral = true => 30-50%
SCENARIO="geometry_v1"   # detector setup
# BFIELD=2               # magnetic field  [kG]
BFIELD=5                 # magnetic field  [kG]

RADIUS=100



# BARRELRAD=20.    # barrel radius      [cm] (right now equal to TOF)
# BARRELLEN=40.    # barrel half length [cm] (right now equal to TOF)
BARRELRAD=100.    # barrel radius      [cm] (right now equal to TOF)
BARRELLEN=280.    # barrel half length [cm] (right now equal to TOF)
# BARRELETA=1.443   # barrel max pseudorapidity
BARRELETA=4.0   # barrel max pseudorapidity
TAILLX=1.0        # tail on left    [q]
TAILRX=1.3        # tail on right   [q]
# set inner TOF parameters
# TOFRAD=20.       # TOF radius      [cm]
# TOFLEN=56.       # TOF half length [cm]
# SIGMAT=0.050  # time resolution [ns]
# # set outer TOF parameters
TOFRAD=100.       # TOF radius      [cm]
TOFLEN=280.       # TOF half length [cm]
SIGMAT=0.020      # time resolution [ns]
SIGMA0=0.200      # vertex time spread [ns]
# set RICH parameters
RICHRAD=100.      # RICH radius      [cm]
RICHLEN=280.      # RICH half length [cm]

### calculate max eta from geometry
# BARRELETA=`awk -v a=$TOFRAD -v b=$TOFLEN 'BEGIN {th=atan2(a,b)*0.5; sth=sin(th); cth=cos(th); print -log(sth/cth)}'`

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


if [[ $SEMICENTRAL = true ]]
  then
    CENTRALITYRANGE="30-50"     # centrality in %
elif [[ $SEMICENTRAL = false ]]
 then
   CENTRALITYRANGE="0-10"       # centrality in %
fi

# copy stuff into running directory
echo " --- selected SYSTEM:     $SYSTEM"

# card
# cp ../delphes/cards/propagate.${BFIELD}kG.tails.tcl propagate.tcl
cp ../delphes/cards/propagate.2kG.tails.tcl propagate.tcl
echo " --- selected B-Field:    ${BFIELD}kG"

# LF and HF weights
cp ./corrWeights/hfe_weights_${CENTRALITYRANGE}central.root hfe_weights.root
cp ./corrWeights/lfe_weights_${CENTRALITYRANGE}central.root lfe_weights.root
echo " --- selected Centrality: ${CENTRALITYRANGE}%"

# resolution files
cp ../resolutionfiles/resolution_test_${BFIELD}kG.root resolution.root
echo " --- selected resolution file: resolution_test_${BFIELD}kG.root"

echo " --- selected SCENARIO:   $SCENARIO"
echo " --- selected RADIUS:     $RADIUS"

# code
cp ./macros/anaEEstudy.h anaEEstudy.h
cp ./macros/anaEEstudy.cxx anaEEstudy.cxx



# LUTS
if [[ $SCENARIO = "werner" ]]
then
  cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.el.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.mu.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.pi.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.ka.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.pr.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.pr.dat
elif [[ $SCENARIO = "default" ]]
then
  cp ../LUTs/lutCovm.${BFIELD}kG.${RADIUS}cm.default/lutCovm.el.${BFIELD}kG.${RADIUS}cm.default.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.${BFIELD}kG.${RADIUS}cm.default/lutCovm.mu.${BFIELD}kG.${RADIUS}cm.default.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.${BFIELD}kG.${RADIUS}cm.default/lutCovm.pi.${BFIELD}kG.${RADIUS}cm.default.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.${BFIELD}kG.${RADIUS}cm.default/lutCovm.ka.${BFIELD}kG.${RADIUS}cm.default.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.${BFIELD}kG.${RADIUS}cm.default/lutCovm.pr.${BFIELD}kG.${RADIUS}cm.default.dat lutCovm.pr.dat
elif [[ $SCENARIO = "geometry_v1" ]]
then
  cp ../LUTs/lutCovm.geometry_v1.rmin${RADIUS}.${BFIELD}kG/lutCovm.el.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.geometry_v1.rmin${RADIUS}.${BFIELD}kG/lutCovm.mu.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.geometry_v1.rmin${RADIUS}.${BFIELD}kG/lutCovm.pi.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.geometry_v1.rmin${RADIUS}.${BFIELD}kG/lutCovm.ka.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.geometry_v1.rmin${RADIUS}.${BFIELD}kG/lutCovm.pr.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.pr.dat
elif [[ $SCENARIO = "its3" ]] # requires a radius of 400cm
then
  cp ../LUTs/lutCovm.its3.rmin${RADIUS}.${BFIELD}kG/lutCovm.el.${BFIELD}kG.rmin${RADIUS}.its3.dat lutCovm.el.dat
  cp ../LUTs/lutCovm.its3.rmin${RADIUS}.${BFIELD}kG/lutCovm.mu.${BFIELD}kG.rmin${RADIUS}.its3.dat lutCovm.mu.dat
  cp ../LUTs/lutCovm.its3.rmin${RADIUS}.${BFIELD}kG/lutCovm.pi.${BFIELD}kG.rmin${RADIUS}.its3.dat lutCovm.pi.dat
  cp ../LUTs/lutCovm.its3.rmin${RADIUS}.${BFIELD}kG/lutCovm.ka.${BFIELD}kG.rmin${RADIUS}.its3.dat lutCovm.ka.dat
  cp ../LUTs/lutCovm.its3.rmin${RADIUS}.${BFIELD}kG/lutCovm.pr.${BFIELD}kG.rmin${RADIUS}.its3.dat lutCovm.pr.dat
else
  echo "!!! check SCENARIO and BFIELD variables"
fi

# make sure B field is set right
sed -i -e "s/set barrel_Bz .*$/set barrel_Bz ${BFIELD}e\-1/" propagate.tcl
sed -i -e "s/double Bz .*$/double Bz = ${BFIELD};/" anaEEstudy.cxx
### set (i)TOF radius
sed -i -e "s/set barrel_Radius .*$/set barrel_Radius ${BARRELRAD}e\-2/" propagate.tcl
sed -i -e "s/double tof_radius = .*$/double tof_radius = ${TOFRAD}\;/" anaEEstudy.cxx
# sed -i -e "s/double inner_tof_radius = .*$/double inner_tof_radius = ${TOFRAD}\;/" anaEEstudy.cxx
### set (i)TOF length
sed -i -e "s/set barrel_HalfLength .*$/set barrel_HalfLength ${BARRELLEN}e\-2/" propagate.tcl
sed -i -e "s/double tof_length = .*$/double tof_length = ${TOFLEN}\;/" anaEEstudy.cxx
# sed -i -e "s/double inner_tof_length = .*$/double inner_tof_length = ${TOFLEN}\;/" anaEEstudy.cxx
### set TOF acceptance
sed -i -e "s/set barrel_Acceptance .*$/set barrel_Acceptance \{ 0.0 + 1.0 * fabs(eta) < ${BARRELETA} \}/" propagate.tcl
### set TOF time resolution and tails
sed -i -e "s/set barrel_TimeResolution .*$/set barrel_TimeResolution ${SIGMAT}e\-9/" propagate.tcl
sed -i -e "s/set barrel_TailRight .*$/set barrel_TailRight ${TAILRX}/" propagate.tcl
sed -i -e "s/set barrel_TailLeft  .*$/set barrel_TailLeft ${TAILLX}/" propagate.tcl
sed -i -e "s/double tof_sigmat = .*$/double tof_sigmat = ${SIGMAT}\;/" anaEEstudy.cxx
sed -i -e "s/double tof_sigma0 = .*$/double tof_sigma0 = ${SIGMA0}\;/" anaEEstudy.cxx
# sed -i -e "s/double inner_tof_sigmat = .*$/double inner_tof_sigmat = ${SIGMAT}\;/" anaEEstudy.cxx
# sed -i -e "s/double inner_tof_sigma0 = .*$/double inner_tof_sigma0 = ${SIGMA0}\;/" anaEEstudy.cxx
### set RICH radius
sed -i -e "s/double rich_radius = .*$/double rich_radius = ${RICHRAD}\;/" anaEEstudy.cxx
### set RICH length
sed -i -e "s/double rich_length = .*$/double rich_length = ${RICHLEN}\;/" anaEEstudy.cxx
### set centrality boolian for anaEEstudy
sed -i -e "s/bool bSemiCentral  = .*$/bool bSemiCentral  = ${SEMICENTRAL}\;/" anaEEstudy.cxx


# adapt pt cuts corresponding to selected B-field
# if [[ $BFIELD -eq 2 ]]
# then
#   PTACC=4 # pt = 40 MeV/c acceptance cut for low B field
# elif [[ $BFIELD -eq 5 ]]
# then
#   PTACC=8 # pt = 40 MeV/c acceptance cut for low B field
# fi
# sed -i -e "s/double PtCut .*$/double PtCut = ${PTACC}e\-2;/" anaEEstudy.cxx
# echo " --- using Pt Acceptance cut:  0.0${PTACC}GeV/c"

#load the right LUTs in the anaEEstudy.cxx

# cp ../preshower/macros/preparePreshowerEff.C .
# root -l -b -q preparePreshowerEff.C
# echo "Finished   preparePreshowerEff.C"
echo ""


echo " ----------------------------------"
echo " List of set detector parameters: "
echo " SIGMAT      = ${SIGMAT}      # time resolution [ns]"
echo " SIGMA0      = ${SIGMA0}      # vertex time spread [ns]"
echo " BARRELRAD   = ${BARRELRAD}       # barrel radius      [cm] (right now equal to TOF)"
echo " BARRELLEN   = ${BARRELLEN}       # barrel half length [cm] (right now equal to TOF)"
echo " BARRELETA   = ${BARRELETA}        # barrel max pseudorapidity"
echo " TAILLX      = ${TAILLX}        # tail on left    [q]"
echo " TAILRX      = ${TAILRX}        # tail on right   [q]"
echo " TOFRAD      = ${TOFRAD}       # TOF radius      [cm]"
echo " TOFLEN      = ${TOFLEN}       # TOF half length [cm]"
echo " RICHRAD     = ${RICHRAD}       # RICH radius      [cm]"
echo " RICHLEN     = ${RICHLEN}       # RICH half length [cm]"
echo " ----------------------------------"
echo ""

DELPHESFILEPATH="./data/delphes_files/$(expr $NEVENTS \* $NRUNS)events_${SYSTEM}_${SCENARIO}_B=${BFIELD}kG"
[ ! -d "$DELPHESFILEPATH" ] &&
(mkdir -p $DELPHESFILEPATH && echo " Created directory: $DELPHESFILEPATH") ||
(echo " Delphes directory already exsits")


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
    (#rm -rf delphes.$I.root &&
     mv delphes.$I.root $DELPHESFILEPATH/delphes.${I}.root
     rm -rf .running.$I && echo " --- run $I completed") ||
    (#rm -rf delphes.$I.root &&
     rm -rf .running.$I && echo " --- run $I crashed") &

done


wait

### Use when one DelphesFile is required
# echo " Merge delphes files..."
# hadd delphes.$(expr $NEVENTS \* $NRUNS).${SYSTEM}.${BFIELD}kG.${CENTRALITYRANGE}central.root delphes.*.root
# mv delphes.$(expr $NEVENTS \* $NRUNS).${SYSTEM}.${BFIELD}kG.${CENTRALITYRANGE}central.root ./data/delphes_files/
# rm -rf delphes.*.root
# echo " Merge delphes files done"





echo " Finished running Delphes "
### merge runs when all done
hadd -f -j $NJOBS anaEEstudy.${SYSTEM}.${SCENARIO}.B=${BFIELD}kG_$(expr $NEVENTS \* $NRUNS)events.root anaEEstudy.*.root &&# rm -rf anaEEstudy.*.root &&
### run analysis scrip on new file
    # cp ./macros/anaPlots.cxx anaPlots.cxx &&
    # root -l -b -q "anaPlots.cxx(\"anaEEstudy.${SYSTEM}.${SCENARIO}.B=0.${BFIELD}_$(expr $NEVENTS \* $NRUNS)events.root\")" &&
    # rm anaPlots.cxx
mv anaEEstudy.${SYSTEM}.${SCENARIO}.B=${BFIELD}kG_$(expr $NEVENTS \* $NRUNS)events.root ./data/prod && rm -rf anaEEstudy.*.root &&

### clean up
rm lutCovm*
rm propagate.tcl
rm delphes.*.log
rm *.root
# rm *.log
rm *.cfg
rm anaEEstudy.cxx
rm anaEEstudy.h
# rm preparePreshowerEff.C
# ./cleanup.sh
