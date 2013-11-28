#!/bin/bash

#### usage: 
# first argument is "hypothesis type". To be chose among: WWELMU,WWELEL,WWMUMU
#
# no 2nd argument --> create cfg.py and submit fwlite macro histoProducer
# 2nd argument="noRun" --> only create cfg.py file
####

EXPECTED_ARGS=2 #at least two arguments
if [ $# -lt $EXPECTED_ARGS ]
then
    echo "you need to provide at least the first argument, which is the hypoType"
    echo "choose between: HWWELMU,HWWMUEL,HWWELEL,HWWMUMU,WWELMU,WWELEL, ..."
    echo "and the pass number: 1,2,3,..."
    exit 1
fi

case $1 in 
    HWWELMU)
	;;
    HWWMUEL)
	;;
    HWWELEL)
	;;
    HWWMUMU)
	;;
    WWELMU)
	;;
    WWELEL)
	;;
    WWMUMU)
	;;
    *)
	echo "ERROR: hypothesis type " $1 "is not known. Exiting"
	exit 1
esac


AnalysisStepFolder=$CMSSW_BASE/src/WWAnalysis/AnalysisStep

IFS=$'\n'; for i in $(cat $AnalysisStepFolder/scripts/listTrees.txt);
do
    if [[ $i == \#* ]] 
    then  continue;
    fi;
    inputStrings=$(echo $i|awk '{print $1}')
    outputString=$(echo $i|awk '{print $2}')
    
    echo "input: " $inputStrings
    echo "output: " $outputString
    label=$(echo $outputString|awk -F. '{print $2}')
    echo "label: " $label

    echo "creating cfg file for sample: " $label;
    sed "s/REPLACE_PASS/00$2/"              $AnalysisStepFolder/test/histoProducer$1_cfg.py > tmp1.py
    sed "s/REPLACEME_INPUT/$inputStrings/"   tmp1.py > tmp2.py
    sed "s/REPLACEME_OUTPUT/$outputString/"  tmp2.py > $label.py; 
    rm tmp1.py;     rm tmp2.py;

    if [[ $2 != "noRun" ]] 
    then
	echo "submitting macro .. "
	histoProducer $label.py >& $label.log &
    fi;
	
done
