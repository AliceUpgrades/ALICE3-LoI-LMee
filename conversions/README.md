### Inputs for the conversion studies.
PYTHIA8 Pb-Pb events are generated with:
```
o2-sim -e TGeant3 -n 5000 -g pythia8hi -m A3IP TRK
```
for cross checks with previous studies we also use hijing with:
```
o2-sim -e TGeant3 -n 5000 -m A3IP TRK -g external --configKeyValues 'GeneratorExternal.fileName=./input/hijing.C'
```
For more details on HIJING
(https://alice.its.cern.ch/jira/browse/AOGM-246)[https://alice.its.cern.ch/jira/browse/AOGM-246]
