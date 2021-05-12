#! /usr/bin/env bash

runDelphes() {
  # write function to be more readable
  ### copy pythia8 configuration and adjust it
  sleep 2 # be save that the file is copied before modifing it. this can get corrupt it when latency is high

  DelphesPythia8 propagate.tcl pythia8.cfg delphes.root &&
  root -q -l "ana.cxx(\"delphes.root\", \"output.root\")"
  #root -q -l "bkg.cxx(\"delphes.root\", \"background.root\")"
}

NEVENTS=$1     # number of events in a run
# SYSTEM="pp_inel"   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.
SYSTEM=$2   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.
ANALYSIS=$3


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
echo " --- generating events:"
echo " --- $NEVENTS $SYSTEM events"

# analyse
cp $ANALYSIS.cxx ana.cxx
# card
cp $DELPHESO2_ROOT/examples/cards/propagate.2kG.tails.tcl propagate.tcl
# code
# pythia configuration
cp pythia8_${SYSTEM}.cfg pythia8.cfg

echo "" >> pythia8.cfg
echo "### run time configuration" >> pythia8.cfg
echo "Main:numberOfEvents $NEVENTS" >> pythia8.cfg
echo "Beams:allowVertexSpread on " >> pythia8.cfg
echo "Beams:sigmaTime 60." >> pythia8.cfg
echo "Random:setSeed on" >> pythia8.cfg
#use process id as seed maybe add time?
echo "Random:seed = $$" >> pythia8.cfg

cp lutCovm.el.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.el.dat
cp lutCovm.mu.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.mu.dat
cp lutCovm.pi.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.pi.dat
cp lutCovm.ka.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.ka.dat
cp lutCovm.pr.werner.rmin${RADIUS}.${BFIELD}kG.dat lutCovm.pr.dat



# Set B fild in propagation card and analysis macro
sed -i -e "s/set barrel_Bz .*$/set barrel_Bz ${BFIELD}e\-1/" propagate.tcl
sed -i -e "s/double Bz .*$/double Bz = ${BFIELD}e\-1;/" ana.cxx

### set TOF radius
sed -i -e "s/set barrel_Radius .*$/set barrel_Radius ${TOFRAD}e\-2/" propagate.tcl
sed -i -e "s/double tof_radius = .*$/double tof_radius = ${TOFRAD}\;/" ana.cxx
### set TOF length
sed -i -e "s/set barrel_HalfLength .*$/set barrel_HalfLength ${TOFLEN}e\-2/" propagate.tcl
sed -i -e "s/double tof_length = .*$/double tof_length = ${TOFLEN}\;/" ana.cxx
### set TOF acceptance
sed -i -e "s/set barrel_Acceptance .*$/set barrel_Acceptance \{ 0.0 + 1.0 * fabs(eta) < ${TOFETA} \}/" propagate.tcl
### set TOF time resolution and tails
sed -i -e "s/set barrel_TimeResolution .*$/set barrel_TimeResolution ${SIGMAT}e\-9/" propagate.tcl
sed -i -e "s/set barrel_TailRight .*$/set barrel_TailRight ${TAILRX}/" propagate.tcl
sed -i -e "s/set barrel_TailLeft  .*$/set barrel_TailLeft ${TAILLX}/" propagate.tcl
sed -i -e "s/double tof_sigmat = .*$/double tof_sigmat = ${SIGMAT}\;/" ana.cxx
sed -i -e "s/double tof_sigma0 = .*$/double tof_sigma0 = ${SIGMA0}\;/" ana.cxx


### create LUTs
# BFIELDT=`awk -v a=$BFIELD 'BEGIN {print a*0.1}'`
# $DELPHESO2_ROOT/examples/scripts/create_luts.sh werner $BFIELDT $TOFRAD



runDelphes &&
(rm -rf delphes.root && echo " --- run $I completed") ||
(rm -rf delphes.root && echo " --- run $I crashed") &

### merge runs when all done
wait

### clean up
# rm lutCovm*
# rm propagate.tcl
# rm *.root
# rm *.log
# rm *.cfg
# rm bkg.cxx
