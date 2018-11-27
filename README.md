# Building the SLntupler

```
cmsrel CMSSW_8_0_28_patch1
git clone https://github.com/NJManganelli/TaggerTest.git Demo
scram b
```

# Running Ntupler
```
cd Demo/ntupler/test
cmsRun fullConfSL.py #Runs all code to produce ntuple, including adding variables (as in below)
cmsRun confSL.py #For case where prepper configuration Demo/prepper/test/D.py has been run instead, producing intermediate output
```
# Running the top tagger
```
cmsrel CMSSW_9_3_12_patch2
git clone git@github.com:susy2015/TopTagger.git
cd TopTagger/TopTagger/test
./configure
make -j8
source taggerSetup.sh
getTaggerCfg.sh -t DeepResolved_DeepCSV_GR_Medium_v1.0.1
#Modify the configuration to lower the discriminant cut to ZERO and change the eta cutoff to 2.4, as training was done without such an eta restriction
cp topTaggerTest.cpp topTaggerTest_original.cpp
ln -s ~/THE/PATH/TO/Demo/ntupler/plugins/.customTopTaggerCode.cpp topTaggerTest.cpp
emacs topTaggerTest.cpp #modify to point to the ntuple created in CMSSW 80X
make
./topTaggerTest
```
