#!/bin/zsh
DIR_INPUT=/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/crab_projects/crab_runTT_test/results
DIR_OUTPUT=crabTT
END_ITER=18
print Input Directory is $DIR_INPUT
print Output Directory is $DIR_OUTPUT
for I in {1..${END_ITER}}; do \.\/topTaggerTest ${DIR_INPUT}\/crabTT_${I}\.root ${DIR_OUTPUT}\/out_${I}\.root > /dev/null; done