#!/bin/zsh
DIR_INPUT=/afs/cern.ch/user/n/nmangane/LTW3/Demo/ntupler/test/crab_projects/crab_runTTSLFilt_test/results
DIR_OUTPUT=crabTTSL
COMMON_FILENAME=crabTTSLFilt
BEGIN_ITER=27
END_ITER=38
print Input Directory is $DIR_INPUT
print Output Directory is $DIR_OUTPUT
for I in {${BEGIN_ITER}..${END_ITER}}; do \.\/topTaggerTest ${DIR_INPUT}\/${COMMON_FILENAME}_${I}\.root ${DIR_OUTPUT}\/o${I}\.root > /dev/null; done