##Full simulations
We use full MC simulations right now to study two topics. The combinatorial background from conversions in the detector material and the tracking efficiency fir the low pt electrons.
The `o2-sim` executable is used in both cases to steer the generation of the MC.
#### Inputs for the conversion studies.
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

#### Low pt tracking
We for now have a collection of scripts to steer the tracking. Thanks to David for sharing his code.
`runbatch.sh`: The "master script" for running this. It needs an argument which is the job number, and it will create an output directory with name job_N (where N is the argument it received). It is like this so that you can conveniently run with stuff like GNU parallel or even some other bash script.

`run_simulation.sh`: The simulation script. You will want to edit this to change generator type and change the number of events being generated and things like  that.

`trigger_multiplicity.C` is an auxiliary file for the generator. We ask that the simulation triggers on the presence of an electron. This can also be used to trigger on the event multiplicity.

`run_strangenesstracking_01.C`: This does the tracking itself and saves a file called "treeouput.root" which saves per-track properties.

`config.ini` and `transport.C`: Nit used right now, should make it possible to select events with certain impact parameters. This could help to implement an easy centrality selection on the generator level.
