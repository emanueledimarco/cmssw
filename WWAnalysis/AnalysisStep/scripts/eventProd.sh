#!/bin/bash


### DOCUMENTATION ###
#
# action:
#  --> submit
######################

action=$1
id=$2

### variable that you want to set ########
DIR="${CMSSW_BASE}/src/WWAnalysis/AnalysisStep"
#pubName=V03_WWEventSkimProd_PU
pubName=V03_WWEventSkimProd_noPU
#pubName=V03_WWEventSkimProd_Nov15_FromSept17ReRecoSkim_v2
#pubName=V03_WWEventSkimProd_Nov15_FromRun2010BPromptReco_v2
#pubName=V03_WWEventSkimProd_FromNov4ReReco_Run2010A
#pubName=V03_WWEventSkimProd_FromNov4ReReco_Run2010B
pset=eventProducer_ALL_cfg
njobs=10
globalTag='START38_V12::All'
##########################################



case $id in

    1)
        ### WZ
	dataset=/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER 
	taskName=id1
	njobs=2
	;;

    2)
        ### ZZ, 244 729 events
	dataset=/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER
	taskName=id2
	njobs=2
	;;

    3)
	### Z+Jets (ee), 630 627 events
	dataset=/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER 
	taskName=id3
	njobs=6
	;;

    4)
	### Z+Jets (mumu), 672 975 events
	dataset=/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER 
	taskName=id4
	njobs=6
	;;

    5)
	### Z+Jets, tautau, 40876 events:
	dataset=/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER
	taskName=id5
	njobs=1
	;;

    6)
	### W+Jets, munu, 47 384
	dataset=/WToMuNu_TuneZ2_7TeV-pythia6/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER 
	taskName=id6
	njobs=1
	;;

    7)
	### W+Jets, elnu, 40047 (partial sample)
	dataset=/WToENu_TuneZ2_7TeV-pythia6/mangano-V03_WWSkim_Nov13_v2-fd9da363ddba42283da31e6938d795fe/USER
	taskName=id7
	njobs=1
	;;

    8)
	### W+Jets, taunu, 10 763
	dataset=/WToTauNu_TuneZ2_7TeV-pythia6-tauola/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER 
	taskName=id8
	njobs=1
	;;

    9)
	### TTbar 2l2un2B, 627 812:
	dataset=/TTTo2L2Nu2B_7TeV-powheg-pythia6/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER 
	taskName=id9
	njobs=4
	;;

    10)
	dataset=/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER
	taskName=id10
	njobs=1
	;;

    11)
	dataset=/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER
	taskName=id11
	njobs=1
	;;

    12)
	dataset=/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER
	taskName=id12
	njobs=1
	;;

    13)
	### WW, 41 910:
	dataset=/WWTo2L2Nu_TuneZ2_7TeV-pythia6/mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER 
	taskName=id13
	njobs=1
	;;

    ### WW, madgraph events:
    14)
	dataset=
	taskName=id14
	njobs=1
	;;

    ### TTbar inclusive events:
    15)
	dataset=
	taskName=id15
	njobs=1
	;;

    ### WW no PU, 110 000 events:
    16)
	dataset=/WWTo2L2Nu_TuneZ2_7TeV-pythia6/mangano-V03_WWSkim_Nov13_noPU-fd9da363ddba42283da31e6938d795fe/USER 
	taskName=id16
	njobs=1
	;;



    #### samples from DATA ###

    61)
	### 
	dataset=/Mu/gpetrucc-Sep17_V03_WWSkim_Nov12-18159a1f6b63e7598ba97d74632060ac/USER
	taskName=id61
	njobs=1
	;;

    62)
	### 
	dataset=/EG/gpetrucc-Sep17_V03_WWSkim_Nov12-18159a1f6b63e7598ba97d74632060ac/USER

	taskName=id62
	njobs=3
	;;

    63)
	### 
	dataset=/Mu/gpetrucc-Run2010B_V03_WWSkim_Nov12-18159a1f6b63e7598ba97d74632060ac/USER
	taskName=id63
	njobs=3
	;;

    64)
	dataset=/Electron/gpetrucc-Run2010B_V03_WWSkim_Nov12-18159a1f6b63e7598ba97d74632060ac/USER
	taskName=id64
	njobs=3
	;;


    65)
	### 
	dataset=
	taskName=id65
	njobs=1
	;;

    66)
	### 
	dataset=

	taskName=id66
	njobs=3
	;;

    67)
	### 
	dataset=
	taskName=id67
	njobs=3
	;;

    68)
	dataset=/Electron/mangano-V03_WWSkim_Nov13-18159a1f6b63e7598ba97d74632060ac/USER
	taskName=id68
	njobs=10
	;;


    *)
	echo -e "WARNING: The id you specified doesn't correspond to any of the datasets explicitely listed\n"
	;;
esac

### Higgs signal datasets
testHiggs=$(($id/100000))
if [[ $testHiggs -ge 1 && $testHiggs -lt 2 ]];
then
    echo -e "GREETINGS: The ID you chose corresponds to a Higgs dataset. ID-Dataset link being determined right now..\n"
    newID=$(($id % 100000)); 
    pu=$(($newID/10000)); echo "pu: " $pu;
    newID=$(($newID % 10000)); 
    sample=$(($newID/1000)); echo "sample: " $sample;
    mass=$(($newID % 1000)); echo "mass: " $mass;
    
    case $pu in
	0)
	    puString=mangano-V03_WWSkim_Nov13-fd9da363ddba42283da31e6938d795fe/USER
	    ;;
	1)
	    puString=
	    ;;
	*)
	    echo "ERROR: pu="$pu
	    echo "For the PU, you have to choose between 0 and 1"
	    ;;
    esac

    case $sample in
	1)
	    sampleString=/GluGluToHToWWTo2L2Nu_M-${mass}_7TeV-powheg-pythia6/
	    njobs=2
	    ;;
	2)
	    sampleString=/GluGluToHToWWToLNuTauNu_M-${mass}_7TeV-powheg-pythia6/
	    njobs=1
	    ;;
	3)
	    sampleString=/GluGluToHToWWTo2Tau2Nu_M-${mass}_7TeV-powheg-pythia6/
	    njobs=1
	    ;;
	4)
	    sampleString=/VBF_HToWWTo2L2Nu_M-${mass}_7TeV-powheg-pythia6/
	    njobs=2
	    ;;
	5)
	    sampleString=/VBF_HToWWToLNuTauNu_M-${mass}_7TeV-powheg-pythia6/
	    njobs=1
	    ;;
	6)
	    sampleString=/VBF_HToWWTo2Tau2Nu_M-${mass}_7TeV-powheg-pythia6/
	    njobs=1
	    ;;
	*)
	    echo "ERROR: sample="$sample
	    echo "For the sample, you have to choose between 1 and 6"
	    ;;
    esac
    # 17 masses: 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600
    dataset=$sampleString$puString
    taskName=id$id

    echo "final dataset name: " $dataset
fi



#

case $action in
    submit) echo "action is submit crab jobs"
	#cat $DIR/test/$pset.py |sed "s#SET_GLOBALTAG#$globalTag#" > pset.$taskName.py       
	cp $DIR/test/$pset.py pset.$taskName.py       

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

