#!/bin/bash


### DOCUMENTATION ###
#
# action:
#  --> submit
#  --> copy
######################

action=$1
id=$2

### variable that you want to set ########
DIR="${CMSSW_BASE}/src/WWAnalysis/SkimStep"
pubName=V03_WWSkim_Nov13_noPU
pset=hwwSkim
njobs=100
#globalTag='START38_V12::All'
##########################################


case $id in
    ### WZ, 2 194 752 events, @T2_UA_KIPT,T2_IT_Rome:
    1)
	dataset=/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO
	taskName=id1
	njobs=150
	;;

    ### ZZ, 21 13 368 events, @T2_UA_KIPT:
    2)
	dataset=/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id2
	njobs=150
	;;


    ### Z+Jets, ee, 1 998 990 events:
    3)
	dataset=/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id3
	njobs=150
	;;

    ### Z+Jets, mm, 1 998 931 events:
    4)
	dataset=/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id4
	njobs=20
	;;

    ### Z+Jets, tautau, ?? events:
    5)
	dataset=/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id5
	njobs=10
	;;

    ### W+Jets, munu, 5 283 540:
    6)
	dataset=/WToMuNu_TuneZ2_7TeV-pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id6
	njobs=400
	;;

    ### W+Jets, enu, 5 021 834:
    7)
	dataset=/WToENu_TuneZ2_7TeV-pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id7-attempt3
	njobs=40
	;;

    ### W+Jets, taunu, 5 193 750:
    8)
	dataset=/WToTauNu_TuneZ2_7TeV-pythia6-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id8
	njobs=400
	;;

    ### TTbar 2l2un2B, 996 022:
    9)
	dataset=/TTTo2L2Nu2B_7TeV-powheg-pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id9
	njobs=80
	;;

    ### SingleTop, s-channel,494 967 events:
    10)
	dataset=/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id10
	njobs=5
	;;

    ### SingleTop, t-channel, 484060 events:
    11)
	dataset=/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id11
	njobs=5
	;;

    ### SingleTop, tW-channel, 494961 events:
    12)
	dataset=/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO 
	taskName=id12
	njobs=5
	;;

    ### WW, 110 000 events:
    13)
	dataset=/WWTo2L2Nu_TuneZ2_7TeV-pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO
	taskName=id13
	njobs=5
	;;

    ### WW, madgraph events:
    14)
	dataset=/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO
	taskName=id14
	njobs=150
	;;

    ### TTbar inclusive events:
    15)
	dataset=/TT_TuneZ2_7TeV-pythia6-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO
	taskName=id15
	njobs=150
	;;

    ### WW no PU, 110 000 events:
    16)
	dataset=/WWTo2L2Nu_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO 
	taskName=id16
	njobs=10
	;;
    
    #### samples from DATA ###

    65)
	dataset=/Mu/Run2010A-Nov4ReReco_v1/RECO
	taskName=id65
	njobs=500
	;;

    66)
	dataset=/EG/Run2010A-Nov4ReReco_v1/RECO
	taskName=id66
	njobs=500
	;;

    67)
	dataset=/Mu/Run2010B-Nov4ReReco_v1/RECO
	taskName=id67
	njobs=250
	;;

    68)
	dataset=/Electron/Run2010B-Nov4ReReco_v1/RECO
	taskName=id68
	njobs=250
	;;


    *)
	echo -e "WARNING: The id you specified doesn't correspond to any of the datasets explicitely listed\n"
	;;
esac



### Higgs signal datasets
testHiggs=$(($id/100000))
if [[ $testHiggs -ge 1 && $testHiggs -lt 2 ]]; 
then  #enter here if id in [100k+1,200k-1]
    echo -e "GREETINGS: The ID you chose corresponds to a Higgs dataset. ID-Dataset link being determined right now..\n"
    newID=$(($id % 100000)); 
    pu=$(($newID/10000)); echo "pu: " $pu;
    newID=$(($newID % 10000)); 
    sample=$(($newID/1000)); echo "sample: " $sample;
    mass=$(($newID % 1000)); echo "mass: " $mass;
    
    case $pu in
	0)
	    puString=Fall10-START38_V12-v1/GEN-SIM-RECO 
	    ;;
	1)
	    puString=Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO
	    ;;
	*)
	    echo "ERROR: pu="$pu
	    echo "For the PU, you have to choose between 0 and 1"
	    ;;
    esac

    case $sample in
	1)
	    sampleString=/GluGluToHToWWTo2L2Nu_M-${mass}_7TeV-powheg-pythia6/
	    ;;
	2)
	    sampleString=/GluGluToHToWWToLNuTauNu_M-${mass}_7TeV-powheg-pythia6/
	    ;;
	3)
	    sampleString=/GluGluToHToWWTo2Tau2Nu_M-${mass}_7TeV-powheg-pythia6/
	    ;;
	4)
	    sampleString=/VBF_HToWWTo2L2Nu_M-${mass}_7TeV-powheg-pythia6/
	    ;;
	5)
	    sampleString=/VBF_HToWWToLNuTauNu_M-${mass}_7TeV-powheg-pythia6/
	    ;;
	6)
	    sampleString=/VBF_HToWWTo2Tau2Nu_M-${mass}_7TeV-powheg-pythia6/
	    ;;
	*)
	    echo "ERROR: sample="$sample
	    echo "For the sample, you have to choose between 1 and 6"
	    ;;
    esac
    # 17 masses: 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600
    dataset=$sampleString$puString
    taskName=id$id
    njobs=10

    echo "final dataset name: " $dataset
fi




case $action in
    submit) echo "action is submit crab jobs"
	#cat $DIR/test/$pset.py |sed "s#SET_GLOBALTAG#$globalTag#" | \
	#    sed "s#SET_OUTPUT#$outputFile#" > pset.$taskName.py
	
	cp $DIR/test/$pset.py  pset.$taskName.py

	cat $DIR/crab/crab.template.cfg |sed "s#SET_DATASET#$dataset#" | \
	    sed "s#SET_PUB_NAME#$pubName#" |sed "s#SET_TASK_NAME#$taskName#" | sed "s#SET_NJOBS#$njobs#" > tmp.cfg
	cat tmp.cfg| sed "s#SET_PSET#$PWD/pset.$taskName.py#" > crab.$taskName.cfg; rm tmp.cfg
	crab -create -cfg crab.$taskName.cfg
	rm crab.$taskName.cfg;
	#rm pset.$taskName.py*;
	crab -submit -c $taskName
	;;

    copy) 
	echo "action = copy is not implemented yet"
	;;


    *) echo "ERROR: input type of action is not recognized. Choose between \"submit\" and \"copy\" "
	exit 1
	;;
esac

