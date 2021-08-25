#! /usr/bin/env bash
# $1, $2 ... are the arguments from the r5lmee_FE.jdl file (Arguments = "$1 $2 $3 $4";)

runDelphes() {
  # write function to be more readable
  ### copy pythia8 configuration and adjust it
  sleep 2 # be save that the file is copied before modifing it. this can get corrupt it when latency is high

  DelphesPythia8 propagate.tcl pythia8.cfg delphes.root &&
  root -q -l "anaEEstudy.cxx(\"delphes.root\", \"anaEEstudy.root\")"
}

NEVENTS=$1     # number of events in a run

# SYSTEM="pp_inel"   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.
SYSTEM=$2   # Select the system. This will copy the coresponding pythia configuration. Make sure it exists in the pythia directory.

RADIUS=100     # radius tracks have to reach for reco
# RADIUS=10     # radius tracks have to reach for reco

BFIELD=$3      # magnetic field  [kG]
SIGMAT=0.020   # time resolution [ns]
# SIGMAT=0.050   # time resolution [ns]
SIGMA0=0.200         # vertex time spread [ns]
TAILLX=1.0     # tail on left    [q]
TAILRX=1.3     # tail on right   [q]
# TOFRAD=19.     # TOF radius      [cm]
# TOFLEN=38.     # TOF half length [cm]
TOFRAD=100.     # TOF radius      [cm]
TOFLEN=200.     # TOF half length [cm]
TOFETA=1.443   # TOF max pseudorapidity
RICHRAD=100.      # RICH radius      [cm]
RICHLEN=200.      # RICH half length [cm]

ENERGY=$4

### calculate max eta from geometry
TOFETA=`awk -v a=$TOFRAD -v b=$TOFLEN 'BEGIN {th=atan2(a,b)*0.5; sth=sin(th); cth=cos(th); print -log(sth/cth)}'`
echo "maxEta = $TOFETA"

#how many events are generated
echo " --- generating events:"
echo " --- $NEVENTS $SYSTEM events"


# card
cp $DELPHESO2_ROOT/examples/cards/propagate.2kG.tails.tcl propagate.tcl
# code
# pythia configuration
cp pythia8_${SYSTEM}_${ENERGY}TeV.cfg pythia8.cfg
# resolution
cp resolution_test_${BFIELD}kG.root resolution.root

echo "" >> pythia8.cfg
echo "### run time configuration" >> pythia8.cfg
echo "Main:numberOfEvents $NEVENTS" >> pythia8.cfg
# echo "Beams:allowVertexSpread on " >> pythia8.cfg
# echo "Beams:sigmaTime 60." >> pythia8.cfg
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
sed -i -e "s/double Bz .*$/double Bz = ${BFIELD};/" anaEEstudy.cxx

### set TOF radius
sed -i -e "s/set barrel_Radius .*$/set barrel_Radius ${TOFRAD}e\-2/" propagate.tcl
sed -i -e "s/double tof_radius = .*$/double tof_radius = ${TOFRAD}\;/" anaEEstudy.cxx
### set TOF length
sed -i -e "s/set barrel_HalfLength .*$/set barrel_HalfLength ${TOFLEN}e\-2/" propagate.tcl
sed -i -e "s/double tof_length = .*$/double tof_length = ${TOFLEN}\;/" anaEEstudy.cxx
### set TOF acceptance
sed -i -e "s/set barrel_Acceptance .*$/set barrel_Acceptance \{ 0.0 + 1.0 * fabs(eta) < ${TOFETA} \}/" propagate.tcl
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


# # adapt pt cuts corresponding to selected B-field
# if [[ $BFIELD -eq 2 ]]
# then
#   PTACC=4 # pt = 40 MeV/c acceptance cut for low B field
# elif [[ $BFIELD -eq 5 ]]
# then
#   PTACC=8 # pt = 40 MeV/c acceptance cut for low B field
# fi
# sed -i -e "s/double PtCut .*$/double PtCut = ${PTACC}e\-2;/" anaEEstudy.cxx
# echo " --- using Pt Acceptance cut:  0.0${PTACC}GeV/c"

### create LUTs
# BFIELDT=`awk -v a=$BFIELD 'BEGIN {print a*0.1}'`
# $DELPHESO2_ROOT/examples/scripts/create_luts.sh werner $BFIELDT $TOFRAD


echo ""
echo " ----------------------------------"
echo " List of set detector parameters: "
echo " SIGMAT      = ${SIGMAT}      # time resolution [ns]"
echo " SIGMA0      = ${SIGMA0}      # vertex time spread [ns]"
echo " BARRELRAD   = ${TOFRAD}       # barrel radius      [cm] (right now equal to TOF)"
echo " BARRELLEN   = ${TOFLEN}       # barrel half length [cm] (right now equal to TOF)"
echo " BARRELETA   = ${TOFETA}    # barrel max pseudorapidity"
echo " TAILLX      = ${TAILLX}        # tail on left    [q]"
echo " TAILRX      = ${TAILRX}        # tail on right   [q]"
echo " TOFRAD      = ${TOFRAD}       # TOF radius      [cm]"
echo " TOFLEN      = ${TOFLEN}       # TOF half length [cm]"
echo " RICHRAD     = ${RICHRAD}       # RICH radius      [cm]"
echo " RICHLEN     = ${RICHLEN}       # RICH half length [cm]"
echo " ----------------------------------"
echo ""


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
# rm anaEEstudy.cxx
