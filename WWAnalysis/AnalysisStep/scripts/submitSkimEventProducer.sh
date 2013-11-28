#!/bin/bash

#### usage: 
# first argument is "hypothesis type". To be chose among: WWELMU,WWELEL,WWMUMU
#
# no 2nd argument --> create cfg.py and submit cmsRun
# 2nd argument="noRun" --> only create cfg.py file
####

EXPECTED_ARGS=2 #at least one argument
if [ $# -lt $EXPECTED_ARGS ]
then
    echo "you need to provide at least the first 2 arguments, which is the hypoType"
    echo "choose between: WWELMU,WWMUEL,WWELEL,WWMUMU"
    echo "and the pass number: 1,2,3,..."
    exit 1
fi

filterEl=False
filterMu=False

case $1 in 
    WWELMU)
	filterEl=True
	filterMu=True
	;;
    WWMUEL)
	filterEl=True
	filterMu=True
	;;
    WWELEL)
	filterEl=True
	filterMu=False
	;;
    WWMUMU)
	filterEl=False
	filterMu=True
	;;
    *)
	echo "ERROR: hypothesis type " $1 "is not known. Exiting"
	exit 1
esac


AnalysisStepFolder=$CMSSW_BASE/src/WWAnalysis/AnalysisStep

IFS=$'\n'
for i in $(cat $AnalysisStepFolder/scripts/listDatasets.txt);
do
    if [[ $i == \#* ]] 
    then  continue;
    fi;
    echo "creating cfg file for sample" $i;
    sed "s/REPLACEME/$i/"  $AnalysisStepFolder/test/eventProducer_cfg.py > tmp1.py
    sed "s/REPLACE_HYPOTYPE/$1/"          tmp1.py > tmp2.py
    sed "s/REPLACE_PASS/$2/"              tmp2.py > tmp3.py
    sed "s/REPLACE_FILTER_EL/$filterEl/"  tmp3.py > tmp4.py
    sed "s/REPLACE_FILTER_MU/$filterMu/"  tmp4.py > $i.py
    rm tmp1.py; rm tmp2.py; rm tmp3.py; rm tmp4.py



    if [[ $2 != "noRun" ]] 
    then
	echo "submitting cmsRun job.. "
	cmsRun $i.py >& $i.log &
    fi;
	
done
