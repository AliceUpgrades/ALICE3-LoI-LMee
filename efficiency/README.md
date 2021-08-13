The sotred histograms to calculate the efficiency can be found in the "rootfiles" stored in this cernbox link:
https://cernbox.cern.ch/index.php/s/XfA79SZaZc5OFKm


To run on the grid use the following files:
anaEEstudy.cxx
generateDelphes_FE.sh
r5lmee_FE.jdl
validation_FE.sh

further make sure that the "resolution files", "LUTs", "corrWeights files" (for LF and HF) and "pythia configs" are uploaded and accessible for the r5lmee.jdl

To submit the Job to the grid use the folling pattern:
NEvents,  coll.System,  BField (kG),  coll.Energy, directoryNameToStoreFiles

An ecample could look like this: 1000 PbPb 5 502 OutputDir
