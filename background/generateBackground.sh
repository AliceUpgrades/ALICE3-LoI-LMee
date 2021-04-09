#! /usr/bin/env bash

runDelphes() {
  ### copy pythia8 configuration and adjust it
  cp ./pythia8.cfg pythia8.$1.cfg

  # write function to be more readable
  echo "Main:numberOfEvents $2" >> pythia8.$1.cfg
  random=$(date +%d%H%M%S)
  # echo "Random:seed = $random" >> pythia8.$1.cfg
  echo "Random:seed = $1" >> pythia8.$1.cfg
  # echo "Beams:allowVertexSpread on " >> pythia8.$I.cfg
  # echo "Beams:sigmaTime 60." >> pythia8.$I.cfg

  DelphesPythia8 propagate.tcl pythia8.$1.cfg delphes.$1.root  &> delphes.$1.log &&
  root -q -l "bkg.cxx(\"delphes.$1.root\", \"background.$1.root\")" &> bkg.$1.log
}

# SYSTEM="pp_inel"   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.
SYSTEM="PbPb"   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.

NJOBS=2        # number of max parallel runs
NRUNS=2     # number of runs

NEVENTS=250    # number of events in a run

RADIUS=100     # radius tracks have to reach for reco

BFIELD=5      # magnetic field  [kG]
SIGMAT=0.020   # time resolution [ns]
SIGMA0=0.200         # vertex time spread [ns]
TAILLX=1.0     # tail on left    [q]
TAILRX=1.3     # tail on right   [q]
TOFRAD=100.    # TOF radius      [cm]
TOFLEN=200.    # TOF half length [cm]
TOFETA=1.443   # TOF max pseudorapidity

### calculate max eta from geometry
TOFETA=`awk -v a=$TOFRAD -v b=$TOFLEN 'BEGIN {th=atan2(a,b)*0.5; sth=sin(th); cth=cos(th); print -log(sth/cth)}'`
echo "maxEta = $TOFETA"

#how many events are generated
ALLEVENTS=$(expr $NEVENTS \* $NRUNS)
echo " --- generating events:"
echo " --- $ALLEVENTS $SYSTEM events"


# card
cp ../delphes/cards/propagate.2kG.tails.tcl propagate.tcl
# code
cp ./macros/bkg.cxx bkg.cxx
# pythia configuration
cp ../pythia/pythia8_${SYSTEM}.cfg pythia8.cfg
# cp $DELPHESO2_ROOT/examples/pythia8/pythia8_inel.cfg pythia8.cfg
# LUTs
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.el.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.el.dat
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.mu.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.mu.dat
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.pi.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.pi.dat
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.ka.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.ka.dat
cp ../LUTs/lutCovm.werner.rmin${RADIUS}.${BFIELD}kG/lutCovm.pr.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.pr.dat

# cp ../../LUTs/default/lutCovm.el.5kG.dat lutCovm.el.dat # default lUTS for test
# cp ../../LUTs/default/lutCovm.mu.5kG.dat lutCovm.mu.dat # default lUTS for test
# cp ../../LUTs/default/lutCovm.pi.5kG.dat lutCovm.pi.dat # default lUTS for test
# cp ../../LUTs/default/lutCovm.ka.5kG.dat lutCovm.ka.dat # default lUTS for test
# cp ../../LUTs/default/lutCovm.pr.5kG.dat lutCovm.pr.dat # default lUTS for test



# Set B fild in propagation card and analysis macro
sed -i -e "s/set Bz .*$/set Bz ${BFIELD}e\-1/" propagate.tcl
sed -i -e "s/double Bz .*$/double Bz = ${BFIELD}e\-1;/" bkg.cxx

### set TOF radius
sed -i -e "s/set barrel_Radius .*$/set barrel_Radius ${TOFRAD}e\-2/" propagate.tcl
sed -i -e "s/double tof_radius = .*$/double tof_radius = ${TOFRAD}\;/" bkg.cxx
### set TOF length
sed -i -e "s/set barrel_HalfLength .*$/set barrel_HalfLength ${TOFLEN}e\-2/" propagate.tcl
sed -i -e "s/double tof_length = .*$/double tof_length = ${TOFLEN}\;/" bkg.cxx
### set TOF acceptance
sed -i -e "s/set barrel_Acceptance .*$/set barrel_Acceptance \{ 0.0 + 1.0 * fabs(eta) < ${TOFETA} \}/" propagate.tcl
### set TOF time resolution and tails
sed -i -e "s/set barrel_TimeResolution .*$/set barrel_TimeResolution ${SIGMAT}e\-9/" propagate.tcl
sed -i -e "s/set barrel_TailRight .*$/set barrel_TailRight ${TAILRX}/" propagate.tcl
sed -i -e "s/set barrel_TailLeft  .*$/set barrel_TailLeft ${TAILLX}/" propagate.tcl
sed -i -e "s/double tof_sigmat = .*$/double tof_sigmat = ${SIGMAT}\;/" bkg.cxx
sed -i -e "s/double tof_sigma0 = .*$/double tof_sigma0 = ${SIGMA0}\;/" bkg.cxx

### create LUTs
# BFIELDT=`awk -v a=$BFIELD 'BEGIN {print a*0.1}'`
# $DELPHESO2_ROOT/examples/scripts/create_luts.sh werner $BFIELDT $TOFRAD

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
hadd -f background.rmin${RADIUS}.${BFIELD}kG.${SYSTEM}.root background.*.root &&
mv background.rmin${RADIUS}.${BFIELD}kG.${SYSTEM}.root ./output/
rm -rf background.*.root &&

### clean up
rm lutCovm*
rm propagate.tcl
rm *.root
# rm *.log
rm *.cfg
rm bkg.cxx
