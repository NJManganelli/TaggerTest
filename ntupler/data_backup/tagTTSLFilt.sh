#!/bin/zsh
DIR_INPUT=/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/crab_projects/crab_runTTSLFilt_test/results
DIR_OUTPUT=crabTTSL
END_ITER=38
print Input Directory is $DIR_INPUT
print Output Directory is $DIR_OUTPUT
for I in {1..${END_ITER}}; do \.\/topTaggerTest ${DIR_INPUT}\/crabTT_${I}\.root ${DIR_OUTPUT}\/o${I}\.root > /dev/null; done