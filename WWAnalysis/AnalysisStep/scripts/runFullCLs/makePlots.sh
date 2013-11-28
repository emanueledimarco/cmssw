#!/usr/bin/env bash

prefix=./datacards.4/hww-4.63fb
nParrallel=7
lumi=4.6
niceLumi=4.6fb
#label=cutBased # or maybe S1 for scenario 1?
label=final # or maybe S1 for scenario 1?


rm -rf plots out 
mkdir -p out plots 




names="comb_0j1j2j comb_2j of_0j sf_0j of_1j sf_1j"
for name in $names; do
	cat limitSummary_final_$name.txt|grep -v -|grep -v obs > tmp.txt
	root -b -l -q PlotLimit.C+"(\"tmp.txt\",\"limit_${niceLumi}_${label}_$name\",\"Hww, L = $lumi fb^{-1}\",   true,true)"
	rm tmp.txt
done





