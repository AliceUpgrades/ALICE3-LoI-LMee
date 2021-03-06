### Background
Estimate the background to the dielectron measurement. This includes for now only the combinatorial background, which is calculated in Pythia8 Angantyr PbPb events.
By running the `generateBackground.sh` script the events are generated and processed with DelphesO2. The outputs are then again processed with `./macros/bkg.C`.
