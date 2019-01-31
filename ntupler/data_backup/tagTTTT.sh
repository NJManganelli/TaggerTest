#!/bin/zsh
DIR_INPUT=/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/crab_projects/crab_runTTTT_test/results
DIR_OUTPUT=crabTTTT
END_ITER=12
print Input Directory is $DIR_INPUT
print Output Directory is $DIR_OUTPUT
for I in {1..${END_ITER}}; do \.\/topTaggerTest ${DIR_INPUT}\/crabTT_${I}\.root ${DIR_OUTPUT}\/o${I}\.root > /dev/null; done