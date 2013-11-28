#!/bin/bash
n=$1
# total events analyzed
cat CMSSW_$n.stdout |grep TrigRe|grep selectionHLT|head -n1 |awk '{print "total : " $4}';

# passing trigger
cat CMSSW_$n.stdout |grep TrigRe|grep selectionHLT|head -n1 |awk '{print "hlt : " $5}';

# passing different cuts
cat CMSSW_$n.stdout |grep TrigRe|grep selection1PairsCounter|head -n1 |awk '{print "cut 1: " $5}'; #1
cat CMSSW_$n.stdout |grep TrigRe|grep ipSelPairsCounter|head -n1 |awk '{print "cut 2: " $5}';  #2
cat CMSSW_$n.stdout |grep TrigRe|grep isoSelPairsCounter|head -n1 |awk '{print "cut 3: " $5}'; #3
cat CMSSW_$n.stdout |grep TrigRe|grep idSelPairsCounter|head -n1 |awk '{print "cut 4: " $5}';  #4
cat CMSSW_$n.stdout |grep TrigRe|grep convSelPairsCounter|head -n1 |awk '{print "cut 5: " $5}'; #5
cat CMSSW_$n.stdout |grep TrigRe|grep metSelLept|head -n1 |awk '{print "cut 6: " $5}';     #6
cat CMSSW_$n.stdout |grep TrigRe|grep diLeptFilter1|head -n1 |awk '{print "cut 7: " $5}';#7
cat CMSSW_$n.stdout |grep TrigRe|grep diLeptSel2|head -n1 |awk '{print "cut 8: " $5}';   #8
cat CMSSW_$n.stdout |grep TrigRe|grep metFilter|head -n1 |awk '{print "cut 9: " $5}';   #9
cat CMSSW_$n.stdout |grep TrigRe|grep jetVetoFilter|head -n1 |awk '{print "cut 10: " $5}';#10
cat CMSSW_$n.stdout |grep TrigRe|grep leptSoftMuonVeto|head -n1 |awk '{print "cut 11: " $5}';#11
cat CMSSW_$n.stdout |grep TrigRe|grep extraLeptonVetoForLept|head -n1 |awk '{print "cut 12: " $5}';#12
