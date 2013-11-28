#/usr/bin/env bash
lumi="4.63"
lumiA="2.12"
lumiB="2.51"
niceLumi="463"
baseDir=$PWD
#masses="`seq 120 10 200` `seq 250 50 600`"
#masses="`seq 120 10 200` `seq 250 50 350`"
#masses="`seq 400 50 600`"
#masses="`seq 120 10 200` `seq 250 50 600`"
masses="110 115 `seq 120 10 200` `seq 250 50 600`"
#masses="110 115"
data=/nfs/bluearc/group/hww/S3/veryCurrentTightTight
#data=/nfs/bluearc/group/hww/S3/currentTightTight
# data=/nfs/bluearc/group/hww/S3/withEfficiencySFs/R42X_S1_V09_S2_V06_S3_V15/TightTight
# data=/nfs/bluearc/group/hww/S3/R42X_S1_V09_S2_V06_S3_V15/TightTight/
juergSys=/afs/cern.ch/user/j/jueugste/public/expUncertainties
s3="$CMSSW_BASE/src/WWAnalysis/AnalysisStep/test/step3/"
dc="$CMSSW_BASE/src/WWAnalysis/AnalysisStep/test/datacard/"
mca="$CMSSW_BASE/src/WWAnalysis/AnalysisStep/python/mcAnalysis.py"
te="$s3/topEstimate.py"
dye="$s3/dyEstimate.py"
nb="$s3/noteBreakdown.py"
c0="$s3/cuts/0JetCuts.txt"
c1="$s3/cuts/1JetCuts.txt"
c2="$s3/cuts/2JetCuts.txt"
eval `for x in $masses; do echo "h$x=$s3/cuts/higgsLevel/$x.txt"; done`
eval `for x in $masses; do echo "b$x=$s3/cuts/bdtLevel/$x.txt"; done`
wwTag=`showtags | grep WWAnalysis | head -n 1 | awk '{print $1}'`

mkdir -p inputCards.$niceLumi && cd inputCards.$niceLumi || die
for M in $masses; do cat $s3/cuts/0JetCuts.txt $s3/cuts/higgsLevel/$M.txt > $M.0j.txt ; done
for M in $masses; do cat $s3/cuts/1JetCuts.txt $s3/cuts/higgsLevel/$M.txt > $M.1j.txt ; done
for M in $masses; do cat $s3/cuts/2JetCuts.txt $s3/cuts/higgsLevel/vbf$M.txt > $M.2j.txt; done
for M in $masses; do 
    python $s3/makeShapeCards.py $s3/mc4mva.txt -2 $M.0j.txt $M.1j.txt $M.2j.txt -P $data --wjAddSyst=0.36 -w --lumi=$lumi --name=cb --cutbased -m $M --dataDr all --bgFolder $dc/final2011 > $M.out 2>&1 &
    pids="$pids $! "
    while [ `ps | grep python | wc -l` -ge 12 ] ; do sleep 10; done
done
wait $pids
# new more better way:
for M in $masses; do
    prefix=hww-${lumi}fb.mH$M
    combineCards.py of0j=$prefix.of_0j_cb.txt sf0j=$prefix.sf_0j_cb.txt > $prefix.comb_0j_cb.txt
    combineCards.py of1j=$prefix.of_1j_cb.txt sf1j=$prefix.sf_1j_cb.txt > $prefix.comb_1j_cb.txt
    combineCards.py of0j=$prefix.of_0j_cb.txt sf0j=$prefix.sf_0j_cb.txt of1j=$prefix.of_1j_cb.txt sf1j=$prefix.sf_1j_cb.txt 2j=$prefix.comb_2j_cb.txt > $prefix.comb_0j1j2j_cb.txt
done

# old way just in case
# for M in $masses; do 
#     #combineCards.py hww-${lumi}fb-cuts.mH$M.elel0j.txt hww-${lumi}fb-cuts.mH$M.mumu0j.txt > hww-${lumi}fb-cuts.mH$M.sf0j.txt
#     #combineCards.py hww-${lumi}fb-cuts.mH$M.elmu0j.txt hww-${lumi}fb-cuts.mH$M.muel0j.txt > hww-${lumi}fb-cuts.mH$M.of0j.txt
#     #combineCards.py hww-${lumi}fb-cuts.mH$M.elel1j.txt hww-${lumi}fb-cuts.mH$M.mumu1j.txt > hww-${lumi}fb-cuts.mH$M.sf1j.txt
#     #combineCards.py hww-${lumi}fb-cuts.mH$M.elmu1j.txt hww-${lumi}fb-cuts.mH$M.muel1j.txt > hww-${lumi}fb-cuts.mH$M.of1j.txt
#     combineCards.py hww-${lumi}fb-cuts.mH$M.[so]f0j.txt > hww-${lumi}fb-cuts.mH$M.comb0j.txt
#     combineCards.py hww-${lumi}fb-cuts.mH$M.[so]f1j.txt > hww-${lumi}fb-cuts.mH$M.comb1j.txt
#     combineCards.py hww-${lumi}fb-cuts.mH$M.comb[012]j.txt > hww-${lumi}fb-cuts.mH$M.comb.txt
# done
cd ..


