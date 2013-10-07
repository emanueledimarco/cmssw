#!/bin/bash

TMP=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH
tar xzf txtFiles.tar.gz
mkdir $CMSSW_BASE/src/Higgs
mkdir $CMSSW_BASE/src/Higgs/Higgs_CS_and_Width
mv txtFiles/ $CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/
tar xzf CreateDatacards.tar.gz
mkdir $CMSSW_BASE/src/HZZ4L_Combination
mkdir $CMSSW_BASE/src/HZZ4L_Combination/CombinationPy
mv CreateDatacards/ $CMSSW_BASE/src/HZZ4L_Combination/CombinationPy/
tar xzf Cards.tar.gz
tar xzf PDFTables.tar.gz
mkdir $CMSSW_BASE/src/ZZMatrixElement/
mkdir $CMSSW_BASE/src/ZZMatrixElement/MEKD/
mkdir $CMSSW_BASE/src/ZZMatrixElement/MEKD/src
mv Cards/ $CMSSW_BASE/src/ZZMatrixElement/MEKD/src/
mv PDFTables/ $CMSSW_BASE/src/ZZMatrixElement/MEKD/src/
cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml -p pset.py
export LD_LIBRARY_PATH=$TMP
