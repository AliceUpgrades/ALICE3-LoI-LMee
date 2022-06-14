#! /usr/bin/env bash

runDelphes() {
  # write function to be more readable
  ### copy pythia8 configuration and adjust it
  cp ./pythia8.cfg pythia8.$1.cfg
  sleep 2 # be save that the file is copied before modifing it. this can get corrupt it when latency is high
  echo "Random:seed = $1" >> pythia8.$1.cfg

  DelphesPythia8 propagate.tcl pythia8.$1.cfg delphes.$1.root  &> delphes.$1.log &&
  root -q -l "doubleTOF.cxx(\"delphes.$1.root\", \"doubleTOF.$1.root\")" &> doubleTOF.$1.log
}

SYSTEM="pp_cc_cr2"   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.
# SYSTEM="PbPb"   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.

NJOBS=25        # number of max parallel runs
NRUNS=25       # number of runs

NEVENTS=10000     # number of events in a run

RADIUS=20     # radius tracks have to reach for reco

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
cp ../../delphes/cards/propagate.2kG.tails.tcl propagate.tcl
# code
cp ./macros/doubleTOF.cxx doubleTOF.cxx
# pythia configuration
cp ../../pythia/pythia8_${SYSTEM}.cfg pythia8.cfg

echo "" >> pythia8.cfg
echo "### run time configuration" >> pythia8.cfg
echo "Main:numberOfEvents $NEVENTS" >> pythia8.cfg
echo "Beams:allowVertexSpread on " >> pythia8.cfg
echo "Beams:sigmaTime 60." >> pythia8.cfg
echo "Random:setSeed on" >> pythia8.cfg



# PY8CFG="pythia8_PbPb"  # pythia8 configuration
# cp $DELPHESO2_ROOT/examples/pythia8/$PY8CFG.cfg pythia8.cfg
# cp $DELPHESO2_ROOT/examples/pythia8/pythia8_inel.cfg pythia8.cfg
# LUTs
cp ../../LUTs/B5kG/lutCovm.el.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.el.dat
cp ../../LUTs/B5kG/lutCovm.mu.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.mu.dat
cp ../../LUTs/B5kG/lutCovm.pi.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.pi.dat
cp ../../LUTs/B5kG/lutCovm.ka.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.ka.dat
cp ../../LUTs/B5kG/lutCovm.pr.${BFIELD}kG.rmin${RADIUS}.geometry_v1.dat lutCovm.pr.dat

# cp ../../LUTs/default/lutCovm.el.5kG.dat lutCovm.el.dat # default lUTS for test
# cp ../../LUTs/default/lutCovm.mu.5kG.dat lutCovm.mu.dat # default lUTS for test
# cp ../../LUTs/default/lutCovm.pi.5kG.dat lutCovm.pi.dat # default lUTS for test
# cp ../../LUTs/default/lutCovm.ka.5kG.dat lutCovm.ka.dat # default lUTS for test
# cp ../../LUTs/default/lutCovm.pr.5kG.dat lutCovm.pr.dat # default lUTS for test



# Set B fild in propagation card and analysis macro
sed -i -e "s/set barrel_Bz .*$/set barrel_Bz ${BFIELD}e\-1/" propagate.tcl
sed -i -e "s/double Bz .*$/double Bz = ${BFIELD}e\-1;/" doubleTOF.cxx

### set TOF radius
sed -i -e "s/set barrel_Radius .*$/set barrel_Radius ${TOFRAD}e\-2/" propagate.tcl
sed -i -e "s/double tof_radius = .*$/double tof_radius = ${TOFRAD}\;/" doubleTOF.cxx
### set TOF length
sed -i -e "s/set barrel_HalfLength .*$/set barrel_HalfLength ${TOFLEN}e\-2/" propagate.tcl
sed -i -e "s/double tof_length = .*$/double tof_length = ${TOFLEN}\;/" doubleTOF.cxx
### set TOF acceptance
sed -i -e "s/set barrel_Acceptance .*$/set barrel_Acceptance \{ 0.0 + 1.0 * fabs(eta) < ${TOFETA} \}/" propagate.tcl
### set TOF time resolution and tails
sed -i -e "s/set barrel_TimeResolution .*$/set barrel_TimeResolution ${SIGMAT}e\-9/" propagate.tcl
sed -i -e "s/set barrel_TailRight .*$/set barrel_TailRight ${TAILRX}/" propagate.tcl
sed -i -e "s/set barrel_TailLeft  .*$/set barrel_TailLeft ${TAILLX}/" propagate.tcl
sed -i -e "s/double tof_sigmat = .*$/double tof_sigmat = ${SIGMAT}\;/" doubleTOF.cxx
sed -i -e "s/double tof_sigma0 = .*$/double tof_sigma0 = ${SIGMA0}\;/" doubleTOF.cxx


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

  runDelphes $I &&
  (rm -rf delphes.$I.root && rm -rf .running.$I && echo " --- run $I completed") ||
  (rm -rf delphes.$I.root && rm -rf .running.$I && echo " --- run $I crashed") &
done

### merge runs when all done
wait
hadd -f doubleTOF.rmin${RADIUS}.${BFIELD}kG.${SYSTEM}.root doubleTOF.*.root &&
mv doubleTOF.rmin${RADIUS}.${BFIELD}kG.${SYSTEM}.root ./output/
rm -rf doubleTOF.*.root &&

### clean up
rm lutCovm*
rm propagate.tcl
rm *.root
#rm *.log
rm *.cfg
rm doubleTOF.cxx
