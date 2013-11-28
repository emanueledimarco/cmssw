#/usr/bin/env bash
lumi="2.12"
niceLumi="212"
baseDir=$PWD
masses="`seq 120 10 200` `seq 250 50 600`"
data=/nfs/bluearc/group/hww/S3/R42X_S1_V07_S2_V04_S3_V10/
eval `for x in $(seq 1 5); do echo "dataS$x=/nfs/bluearc/group/hww/S3/R42X_S1_V07_S2_V05_S3_V11/Scenario.$x/;"; done`
topData=/nfs/bluearc/group/hww/S3/tempTreesForBtagging.R42X_S1_V06_S2_V04_S3_V07++/
s3="$CMSSW_BASE/src/WWAnalysis/AnalysisStep/test/step3/"
mca="$CMSSW_BASE/src/WWAnalysis/AnalysisStep/python/mcAnalysis.py"
te="$s3/topEstimate.py"
dye="$s3/dyEstimate.py"
nb="$s3/noteBreakdown.py"
ema="/afs/cern.ch/user/e/emanuele/www/Analysis/HtoWW/7TeV/LP11BKG"
c0="$s3/cuts/0JetCuts.txt"
c1="$s3/cuts/1JetCuts.txt"
c2="$s3/cuts/2JetCuts.txt"
eval `for x in $masses; do echo "h$x=$s3/cuts/higgsLevel/$x.txt"; done`
eval `for x in $masses; do echo "b$x=$s3/cuts/bdtLevel/$x.txt"; done`

#data card info
svndatacardrepo='svn+ssh://svn.cern.ch/reps/cmshcg/trunk/eps'
guilafsdatacards='/afs/cern.ch/user/c/ceballos/public/limits/inputLimits/'

# Skim Step (S1)
# Make sure to change the VV filter

# Step 2
# Have to worry about triggers

# Step 3 @ UAF:
# Have to worry about jsons

cd $s3
id='m*/*/R42X_S1_V05_S2_V02_S3_V02_ID%(id)s_*/*/*.root'
had=/hadoop/cms/store/user/
# could do it this way but you have to make sure that
# your pattern doesnt pick up anything random (v2 for e.g.)
# python makeStep3.py -p $id allBackgroundSamples $had
# python makeStep3.py -p $id allSignalSamples $had
# python makeStep3.py -p $id -j certifiedLatinos.42X data42xSamples $had
python makeStep3.py allBackgroundSamples 
python makeStep3.py allSignalSamples 
python makeStep3.py -j certifiedLatinos.42X data42xSamples 
bulkcp -v  -r uaf-8.t2.ucsd.edu -n 9 /home/users/mlebourgeois/higgs/hww/master/CMSSW_4_2_7_patch1/src/WWAnalysis/AnalysisStep/test/step3/*.root .

#  data yields
mkdir -p dataYields && cd dataYields || die
for j in 0 1 2 ; do  
    python $nb -w -c 10 -f -T dataYields_${j}j.twiki -L dataYields_${j}j.tex mc4mva.txt $s3/cuts/${j}JetCuts.txt > /dev/null 2>&1 & 
done

# expected yields:
mkdir -p expectedTables && cd expectedTables || die

# convert to the correct headings (Future TODO remove the SFs)
sed -e '/^# *$/d' -e '/^#.*/d' -e 's/ggH     /$m_H=%d$/'  -e 's/vbfH    /$m_H=%d$/' -e 's/VV   \( *\):/WZ+ZZ\1:/'       \
    -e 's/DYEE        /Z$\/\\gamma^*$/'  -e 's/DYMM        /Z$\/\\gamma^*$/' -e 's/DYTT        /Z$\/\\gamma^*$/'           \
    -e '/data/d'               -e 's/^Top     \( *: 01[1-8]\)/single t\1/'   -e 's/^Top       /$t\\bar{t}$/'          \
    -e 's/^WJet  /W+jets/'       $s3/mc4mva.txt > temp.txt


# build the correct order for the WW level summary tables
grep  "03[0-9]_" temp.txt > mc4WWSummary.txt
grep  "01[0-9]_" temp.txt >> mc4WWSummary.txt
grep  "08[0-9]_" temp.txt >> mc4WWSummary.txt
grep  "07[0-9]_" temp.txt >> mc4WWSummary.txt
grep  "001_" temp.txt >> mc4WWSummary.txt
grep  "000_" temp.txt >> mc4WWSummary.txt
grep  "0[0-9][0-9]_" temp.txt | sed 's/^.*\( *\):/tot. bkg.     :/' >> mc4WWSummary.txt
for x in 120 130 140 160 180 200 250 400; do 
    grep  -v "0[0-9][0-9]_" temp.txt | sed "s/%d/$x/g" >> mc4WWSummary.txt
done

# build the correct order for the mass dependent summary tables
grep  -v "m_H" mc4WWSummary.txt > mc4HWWSummary.txt
grep     "m_H"  temp.txt       >> mc4HWWSummary.txt

# done with you
rm temp.txt

for newData in $dataS1 $dataS2 $dataS3 $dataS4; do 
    mkdir -p `basename $newData`
    cd `basename $newData`
    mkdir -p tex twiki
    expectedYields=tex/expectedYields_$niceLumi.tex
    rm $expectedYields
    for j in 0 1 2 ; do
        echo "\\clearpage \\section{${j} Jet Bin}" >> $expectedYields
        caption="Expected yields from simulation for $lumi \invfb at the mass-independent level in the $j-jet bin"
        python $nb -P $newData -w -l $lumi -C "$caption" -c 9 --size="scriptsize" -f -T twiki/wwSummary_${j}j_$niceLumi.twiki  \
            -L tex/wwSummary_${j}j_$niceLumi.tex ../mc4WWSummary.txt $s3/cuts/${j}JetCuts.txt > /dev/null 2>&1 &
        while [ `ps | grep python | wc -l` -ge 7 ] ; do sleep 10; done
        echo "\\input{tex/wwSummary_${j}j_$niceLumi}" >> $expectedYields
    done
    for j in 0 1 2 ; do
        for m in $masses; do 
            caption="Expected yields from simulation for $lumi \invfb for a mass hypothosis of $m GeV in the $j-jet bin"
            [ $j -eq 2 ] && massSpecific="higgsLevel/vbf.txt" || massSpecific="higgsLevel/$m.txt"
            python $nb -P $newData -w -l $lumi -C "$caption" -c 10 --size="scriptsize" -f -m $m -T twiki/hwwSummary_m${m}_${j}j_$niceLumi.twiki  \
                -L tex/hwwSummary_m${m}_${j}j_$niceLumi.tex ../mc4HWWSummary.txt $s3/cuts/${j}JetCuts.txt $s3/cuts/$massSpecific > /dev/null 2>&1 &
            while [ `ps | grep python | wc -l` -ge 7 ] ; do sleep 10; done
            echo "\\input{tex/hwwSummary_m${m}_${j}j_$niceLumi}" >> $expectedYields
        done
    done
    cd ..
done

# update the scale factors:

# Wjets
echo "Starting Wjet Estimate"
mkdir -p wjet && cd wjet || die
echo "No estimate available yet :-(" | tee estimate.$lumi.txt
echo "Zero jet yields: " | tee yields.$lumi.txt
python $mca -w -l $lumi $data/tree_WJetsToLNuMad.root $s3/cuts/0JetCuts.txt | tee -a yields.$lumi.txt
echo "One jet yields: " | tee -a yields.$lumi.txt
python $mca -w -l $lumi $data/tree_WJetsToLNuMad.root $s3/cuts/1JetCuts.txt | tee -a yields.$lumi.txt
grep b-tag yields.$lumi.txt | cut -c 24-31,44-51,64-71,84-91 > sfWorksheet
echo -e "\nYou should replace the last two lines of sfWorksheet with the data yields and then run the following three:
     mcWjet=\`head -n 2 sfWorksheet | awk '{a+=\$1;b+=\$2;c+=\$3;d+=\$4} END {print a,b,c,d}'\`
     dataWjet=\`tail -n 2 sfWorksheet | awk '{a+=\$1;b+=\$2;c+=\$3;d+=\$4} END {print a,b,c,d}'\`
     echo SFs: \$(for x in 1 2 3 4 ; do calc \`echo \$dataWjet | c \$x\`/\`echo \$mcWjet | c \$x\`; done) >> yields.$lumi.txt"
cd $baseDir

# Drell-Yan
echo "Starting Drell-Yan estimates"
mkdir -p dy && cd dy || die
python $dye -f -P $data -l $lumi -w            -O normWW.root  $s3/cuts/0JetCuts.txt > estimate_$lumi.txt        2>&1 &
python $dye -f -P $data -l $lumi -w  --subVV   -O normWW.root  $s3/cuts/0JetCuts.txt > estimate_subVV_$lumi.txt  2>&1 &
cd $baseDir
# Still have to fix this script to be used correctly with the new multiple cuts files
# for x in $s3/cuts/mT/cuts*; do
#     outName=norm`basename $x | grep -o "[0-9][0-9][0-9]"`
#     python $dye --42X -P $data -l $lumi -w --madgraph -O $outName.root $s3/cuts/0JetCuts.txt $x > $outName.txt 2>&1 &
#     while [[ `ps | grep python | grep -v grep | wc -l` -ge 8 ]] ; do sleep 10s; done
#     for y in subVV; do
#         outName=$y`basename $x | grep -o "[0-9][0-9][0-9]"`
#         python $dye --42X -P $data -l $lumi -w --$y --madgraph -O $outName.root $s3/cuts/0JetCuts.txt $x > $outName.txt 2>&1 &
#         while [[ `ps | grep python | grep -v grep | wc -l` -ge 8 ]] ; do sleep 10s; done
#     done
# done

#TTbar
# temp: python topEstimate.py -f -w -l $lumi -P $data mc4mva.txt $s3/cuts/0JetCuts.txt
echo "Starting Top Estimate"
mkdir -p top && cd top || die
python $te -w -P $dataS1 -l $lumi -f -i $s3/mc4mva.txt $c0 | tee estimate_$niceLumi.txt 
set `grep data/mc estimate_$niceLumi.txt `
sf0=$3
err0=$5
sf1=$8
err1=${10}
set `grep '\^{control} ' estimate_$niceLumi.txt `
gam0=$3
gam1=$6
echo $sf0 $err0 $sf1 $err1 $gam0 $gam1
python $s3/quickYields.py $s3/mc4mva.txt Top $c0 $s3/cuts/higgsLevel/%d.txt --jets 0 --scale $sf0 $err0 -w --masses $s3/masses --lumi=$lumi -P $dataS1 --gamma=$gam0
python $s3/quickYields.py $s3/mc4mva.txt Top $c1 $s3/cuts/higgsLevel/%d.txt --jets 1 --scale $sf1 $err1 -w --masses $s3/masses --lumi=$lumi -P $dataS1 --gamma=$gam1
cd $baseDir

# WW
echo "Starting WW Estimate"
mkdir ww && cd ww || die
echo "No estimate available yet :-(" | tee estimate_$lumi.txt
cd $baseDir


# Make some plots:
echo "Have you updated the scale factors?"
mkdir -p scaledPlots && cd scaledPlots || die
sed -e 's/data/Data/'   -e 's/ggH/signal/'  \
    -e 's/vbfH/signal/' -e 's/VV:/WZ:/'       \
    -e 's/DYEE/ZJets/'  -e 's/DYMM/ZJets/'  \
    -e 's/DYTT/ZJets/'  -e '/^#/d'          \
    -e '/Vg/d'          -e 's/WJet:/WJets:/'  \
    -e 's/ggWW/WW/' $s3/mc4mva.txt > $s3/mc4Plots.txt
pids=""
for mass in 130 160 200 ; do 
    python $s3/plotMakerNice.py -f -l $lumi -m $mass -w -P $data --mva=BDTG_boris_$mass -O 0jWW $s3/mc4Plots.txt $s3/cuts/0JetCuts.txt $s3/wwPlots.txt > /dev/null 2>&1 &
    pids="$pids $!"                
    python $s3/plotMakerNice.py -f -l $lumi -m $mass -w -P $data --mva=BDTG_boris_$mass -O 1jWW $s3/mc4Plots.txt $s3/cuts/1JetCuts.txt  $s3/wwPlots.txt > /dev/null 2>&1 &
    pids="$pids $!"
done
wait $pids
cpPlotsToAfs.sh 0jWW hww/plots/$niceLumi/0jWW  &
cpPlotsToAfs.sh 1jWW hww/plots/$niceLumi/1jWW  &
cd ~/afs/www/main/hww/plots/$niceLumi/
for x in `find . -name "*.eps" -print`; do convert $x ${x/eps/pdf}; done
cd $baseDir

# Start here for limit calculation

# cp $ema/*/*.txt $s3/../datacard/
# mv DYCard_mmsplit_0j.txt DYCard_mm_0j.txt
# mv DYCard_mmsplit_1j.txt DYCard_mm_1j.txt
# mv DYCard_eesplit_1j.txt DYCard_ee_1j.txt
# mv DYCard_eesplit_0j.txt DYCard_ee_0j.txt
# perl $s3/../datacard/text2card-WJet.pl 

# Convert them to limit cards
# MANGANO
# input DD cards must be in WWAnalysis/AnalysisStep/test/datacard/
# possibly to change WWAnalysis/AnalysisStep/python/datacardWriter.py
# 
mkdir -p inputCards && cd inputCards || die
for M in $masses; do cat $s3/cuts/0JetCuts.txt $s3/cuts/higgsLevel/$M.txt > $M.txt ; done
for M in $masses; do 
    python $s3/makeShapeCards.py $s3/mc4mva.txt $M.txt  -P $data -w --lumi=$lumi --name=cuts --cutbased -m $M --dataDr all > $M.out 2>&1 &
done

# run the asymptotic limits
# MANGANO
mkdir -p asymptoticLimits && cd asymptoticLimits || die
# edit the variables at the beginning of runTheLimits.sh
source runTheLimits.sh
# output plots in plots

# and run?!?
# 00:55:28 Giovanni: so, if you have on AFS a CMSSW area that contains HiggsAnalysis/CombinedLimit package
# 00:55:52 Giovanni: you can grab from /home/gpetrucc/scratch0/higgs/CMSSW_4_1_4/src/WWAnalysis/AnalysisStep/test/datacard
# 00:56:03 Giovanni: the mk* scripts
# 00:56:54 Giovanni: and from ~gpetrucc/scratch0/higgs/CMSSW_4_1_4/src/HiggsAnalysis/CombinedLimit/test/  the combine*sh
# 00:57:23 Giovanni: fix the paths within them, and then if you run mkUnMilioneDiJobs > jobs; bash jobs
# 00:57:29 Giovanni: it will submit the jobs for the limits
# ... combine /afs/cern.ch/user/m/mwlebour/scratch0/higgs/hww/CMSSW_4_2_4_patch2/src/cards/hww-0.7105fb-cuts.mH120.comb.txt -m 120 -n HWW_comb -M MarkovChainMC -H ProfileLikelihood -t 10 -v -1 -s 202
# hadd2 -f C_154fbema.root 162*/higgsCombineHWW_comb.MarkovChainMC.mH???.*.root
# hadd2 -f C_154fbema_0j.root 162*/higgsCombineHWW_comb_0j.MarkovChainMC.mH???.*.root
# hadd2 -f C_154fbema_1j.root 162*/higgsCombineHWW_comb_1j.MarkovChainMC.mH???.*.root
# 
# cd hww-bogus
# root.exe -b -l -q ../../test/plotting/bandUtils.cxx+ makeBands.cxx
# root.exe -b -l -q ../../test/plotting/bandUtils.cxx+ tdrstyle.C hwwplots.cxx
# 
# # here are the new links for runnign the limit calc
# # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit#HybridNew_algorithm_usage_for_co
# 
# # example for attaching mva:
# python $s3/plotMakerNice.py -f -l 1.091 -m $mass -w -P $data --mva=BDTG_boris_$mass -O 0j mc4EPS.txt 0JetCuts.txt wwPlots.txt 
# 
# # example for attaching the lepW and using the proper reweighting
# $mca  -w -W "lepW*baseW*puW*kfW" --eff=lepW  /nfs/bluearc/group/hww/S3/R42X_S1_V06_S2_V03_S3_V06/tree_031_42X_DYtoMuMu.root $s3/cuts/1JetCuts.txt $s3/cuts/mT/cuts200.txt  -m 200  -f -l 1.54
# 
