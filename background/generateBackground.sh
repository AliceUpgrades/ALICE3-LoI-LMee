#! /usr/bin/env bash

runDelphes() {
  ### copy pythia8 configuration and adjust it
  cp ../pythia/pythia8_PbPb.cfg pythia8.$1.cfg
  # write function to be more readable
  echo "Main:numberOfEvents $2" >> pythia8.$1.cfg
  random=$(date +%d%H%M%S)
  echo "Random:seed = $random" >> pythia8.$1.cfg

  DelphesPythia8 propagate.tcl pythia8.$1.cfg delphes_PbPb.$1.root  &> delphes_PbPb.$1.log &&
  root -q -l "bkg.cxx(\"delphes_PbPb.$1.root\", \"background.$1.root\")" &> bkg.$1.log
}

SYSTEM="PbPb"   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.

NJOBS=1        # number of max parallel runs
NRUNS=1     # number of runs

NEVENTS=100    # number of events in a run

BFIELD=2       # magnetic field  [kG]
RADIUS=100     # radius tracks have to reach for reco

#how many events are generated
ALLEVENTS=$(expr $NEVENTS \* $NRUNS)
echo " --- generating events:"
echo " --- $ALLEVENTS $SYSTEM events"


# card
cp ../delphes/cards/propagate.2kG.tcl propagate.tcl
# code
cp ./macros/bkg.cxx bkg.cxx
# pythia configuration
cp ../pythia/pythia8_${SYSTEM}.cfg pythia8.cfg

# LUTs
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.el.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.el.dat
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.mu.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.mu.dat
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.pi.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.pi.dat
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.ka.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.ka.dat
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.pr.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.pr.dat

# Set B fild in propagation card and analysis macro
sed -i -e "s/set Bz .*$/set Bz ${BFIELD}e\-1/" propagate.tcl
sed -i -e "s/double Bz .*$/double Bz = ${BFIELD}e\-1;/" bkg.cxx

for (( I = 1; I <= $NRUNS; I++ )); do
  while [ $(ls .running.* 2>/dev/null | wc -l) -ge $NJOBS ]; do
    echo " --- waiting for a free slot"
    sleep 10
  done

  ### book the slot
  echo " --- starting run $I"
  touch .running.$I

  runDelphes $I $NEVENTS &&
  (rm -rf delphes.$I.root && rm -rf .running.$I && echo " --- run $I completed") ||
  (rm -rf delphes.$I.root && rm -rf .running.$I && echo " --- run $I crashed") &
done

### merge runs when all done
wait
hadd -f background.rmin${RADIUS}.${BFIELD}kG.root background.*.root &&
mv background.rmin${RADIUS}.${BFIELD}kG.root ./output
rm -rf background.*.root &&

### clean up
rm lutCovm*
rm propagate.tcl
rm *.root
rm *.log
rm *.cfg
rm bkg.cxx
