#!/bin/bash

#for x in ???/0*; do 
#    cd $x
#    echo $x `crab -status | grep ' Total Jobs' | c 2` `ls -1 crab_0_*/res/*.root | wc -l`
#    cd ../.. 
#done > rootFileSummary

cat rootFileSummary | while read x; do
    dir="`echo $x | c 1`/crab_0_*/res/"
    fn=`echo $x | c 1 | awk -F. '{print $2}'`
    nJobs=`echo $x | c 2`
    for i in `seq 1 $nJobs`; do 
        nFiles=`ls -1 $dir/${fn}_${i}_* | wc -l`
        rm `ls -1 $dir/${fn}_${i}_* | head -n $((nFiles-1))`
    done
done
    
