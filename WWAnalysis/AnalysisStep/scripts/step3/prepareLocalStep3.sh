#!/usr/bin/env bash

wwTag=`showtags | grep WWAnalysis | head -n 1 | awk '{print $1}'`
#wwTag=stiCazzi
jsonFile="$CMSSW_BASE/src/WWAnalysis/Misc/Jsons/certifiedUCSD.json"
inputDirectory=/nfs/bluearc/group/trees/hww/R414_S1_V06_S2_V03_Reduced/
mass=${1:-160}
echo $HOSTNAME | grep ucsd && prefix="" || prefix="file:"
    
cat loopOverAll.master | while read x; do
    if [ "${x:0:1}" = "#" ] ; then continue;fi
    echo $x | grep "10[1-6][0-9][0-9][0-9]" > /dev/null 2>&1 
    if [[ $? -eq 0 ]] ; then
        echo $x | grep "10[1-6]$mass" > /dev/null 2>&1 || continue
    fi
    sampleNum=`echo $x | c 1`
    sampleName=`echo $x | c 2`
    numEvents=`echo $x | c 3`
    numJobs=`echo $x | c 4`
    isMC=`echo $x | c 5`
    isGluGlu=`echo $x | c 6`
    sampleType=`echo $x | c 7`
    isLocal="True"

    inputFilesTmp=`ls $inputDirectory/$sampleNum/*.root |awk '{print "\"file:"$1"\","}'`
    inputFiles=`echo $inputFilesTmp|tr -d " "`
    #echo "inputFiles: " $inputFiles


    mkdir -p m$mass/$sampleNum.$sampleName

    cat createEventHistsFullFW_cfg.py | \
        sed s/RMMEFN/$sampleName.root/g |  \
        sed s/RMMEMC/$isMC/g |  \
        sed s/RMMENUM/$sampleNum/g | \
        sed s/RMMENAME/$sampleName/g | \
        sed s/RMMEMASS/$mass/g | \
        sed s/RMMESAMPLE/$sampleType/g | \
        sed s/RMMEGLUGLU/$isGluGlu/g | \
        sed s#RMMEINPUTFILES#$inputFiles#g | \
        sed s/RMMEISLOCAL/$isLocal/g | \
        sed s/RMMEJSON/${jsonFile//\//\\/}/g > m$mass/$sampleNum.$sampleName/$sampleName.py

done
