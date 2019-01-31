#!/bin/zsh
DIR_INPUT=/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/crab_projects/crab_runTT_test/results
DIR_OUTPUT=crabTT
COMMON_FILENAME=crabTT
BEGIN_ITER=1
END_ITER=18
print Input Directory is $DIR_INPUT
print Output Directory is $DIR_OUTPUT
for I in {${BEGIN_ITER}..${END_ITER}}; do \.\/topTaggerTest ${DIR_INPUT}\/${COMMON_FILENAME}_${I}\.root ${DIR_OUTPUT}\/out_${I}\.root > /dev/null; done