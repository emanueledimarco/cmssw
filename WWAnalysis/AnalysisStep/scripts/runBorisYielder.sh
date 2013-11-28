#!/bin/bash


### settings ###
inputFolder=$1
isOnCastor=False


for mode in mumu muel elmu elel 
do
    mode2=""
    case "$mode" in
	mumu)  mode2=MUMU ;;
	muel)  mode2=MUEL ;;
	elmu)  mode2=ELMU ;;
	elel)  mode2=ELEL ;;
    esac
    
    config=$CMSSW_BASE/src/WWAnalysis/AnalysisStep/test/borisYieldProducerHWW${mode2}_cfg.py
    sed "s/label/ww${mode}${label}/" $config | sed "s#replace_inputFolder#${inputFolder}#" | sed "s/replace_isOnCastor/$isOnCastor/" \
    >  config.${mode}.py
    borisYieldProducer config.${mode}.py >& ${mode}.log &
done

