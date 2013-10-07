#!/usr/bin/env bash

prefix=./hww-4.63fb-cuts
nParrallel=7
lumi=4.6
niceLumi=4.6fb
label=cutBased # or maybe S1 for scenario 1?
jets="0j 1j"
#channels="comb elel elmu muel mumu"
#channels="comb elel elmu muel mumu sf of"
channels="sf of"


masses="`seq 120 10 200` `seq 250 50 350`"
mkdir -p out plots root
pids=
for mass in $masses; do 
    combine -M Asymptotic  --rMax 9.0 -v 2 -m $mass -n _${label}_comb $prefix.mH$mass.comb.txt 2> /dev/null | tee out/$mass.comb.out > /dev/null &
    pids="$! $pids"
    #uncomment me when you have the 2j cards
    for j in $jets ; do
	for ch in $channels; do
            combine -M Asymptotic  --rMax 9.0 -v 2 -m $mass -n _${label}_${ch}_$j $prefix.mH$mass.${ch}${j}.txt 2>&1 | tee out/$mass.${ch}_$j.out > /dev/null &
            pids="$! $pids"
            while [ `ps | grep combine | wc -l` -ge $nParrallel ] ; do sleep 10; done
	done
    done
done
wait $pids
mv *.root root/

rm limitSummary_${label}_comb.txt
for mass in $masses; do
    a=`grep "Limit:"                       out/$mass.comb.out | awk '{print $5}'`
    b=`grep "Median for expected limits:"  out/$mass.comb.out | awk '{print $5}'`
    c=`grep "Expected  2.5%:"              out/$mass.comb.out | awk '{print $5}'`
    d=`grep "Expected 16.0%:"              out/$mass.comb.out | awk '{print $5}'`
    e=`grep "Expected 84.0%:"              out/$mass.comb.out | awk '{print $5}'`
    f=`grep "Expected 97.5%:"              out/$mass.comb.out | awk '{print $5}'`
    echo $mass $a 99. $b $c $d $e $f >> limitSummary_${label}_comb.txt
done


for j in $jets; do
    for ch in $channels; do
	rm limitSummary_${label}_${ch}_$j.txt
	for mass in $masses; do
            a=`grep "Limit:"                       out/$mass.${ch}_$j.out | awk '{print $5}'`
            b=`grep "Median for expected limits:"  out/$mass.${ch}_$j.out | awk '{print $5}'`
            c=`grep "Expected  2.5%:"              out/$mass.${ch}_$j.out | awk '{print $5}'`
            d=`grep "Expected 16.0%:"              out/$mass.${ch}_$j.out | awk '{print $5}'`
            e=`grep "Expected 84.0%:"              out/$mass.${ch}_$j.out | awk '{print $5}'`
            f=`grep "Expected 97.5%:"              out/$mass.${ch}_$j.out | awk '{print $5}'`
            echo $mass $a 99. $b $c $d $e $f >> limitSummary_${label}_${ch}_$j.txt
	done;
    done;
done;

                                                                                                          #don't peak!
root -b -l -q PlotLimit.C+"(\"limitSummary_${label}_comb.txt\",   \"limit_${niceLumi}_${label}_comb\",\"Hww, combined, L = $lumi fb^{-1}\",true,true)"
for j in $jets; do
    for ch in $channels; do
	root -b -l -q PlotLimit.C+"(\"limitSummary_${label}_${ch}_$j.txt\",\"limit_${niceLumi}_${label}_${ch}_$j\",\"Hww, 0-jet, L = $lumi fb^{-1}\",   true,true)"
    done
done




