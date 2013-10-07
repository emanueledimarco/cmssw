import FWCore.ParameterSet.Config as cms
# to remove lines with no dataset:
# %s/^[0-9]\+ *$//
# convert rest of the lines to dictionaries
# :%s/^\([0-9]\+\) *\(.*\)$/    '\1':'\2',/

stepOneDatasets = {
    '001':'/GluGluToWWTo4L_TuneZ2star_8TeV-gg2ww-pythia6/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID001_GluGluToWWTo4L-bc58659808d390d23203e20a37543a0e/USER',
	'000':'/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID000_WWJets2LMad-bc58659808d390d23203e20a37543a0e/USER',

	'100':'/SingleElectron/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID100_SingleElectron2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'101':'/SingleMu/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID101_SingleMuon2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'102':'/DoubleElectron/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID102_DoubleElectron2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'103':'/DoubleMu/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID103_DoubleMuon2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'104':'/MuEG/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID104_MuEG2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',

	'200':'/SingleElectron/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID200_SingleElectron2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'201':'/SingleMu/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID201_SingleMuon2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'202':'/DoubleElectron/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID202_DoubleElectron2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'203':'/DoubleMu/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID203_DoubleMuon2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'204':'/MuEG/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID204_MuEG2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',

	'012':'/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID012_TbartWFullDR-bc58659808d390d23203e20a37543a0e/USER',
	'011':'/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID011_TtWFullDR-bc58659808d390d23203e20a37543a0e/USER',
	'010c1':'/TTJets_TuneZ2star_8TeV-madgraph-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID010_TTJetsMadChunck1-55896a93bc94534beab0b43cf4bc1f9b/USER',
        '010c2':'/TTJets_TuneZ2star_8TeV-madgraph-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID010_TTJetsMadChunck2-55896a93bc94534beab0b43cf4bc1f9b/USER',

	'080':'/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID080_WJetsToLNuMad-bc58659808d390d23203e20a37543a0e/USER',

	'070':'/WZ_TuneZ2star_8TeV_pythia6_tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID070_WZFull-bc58659808d390d23203e20a37543a0e/USER',
	'077':'/ZZTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID077_ZZTo2L2Nu-bc58659808d390d23203e20a37543a0e/USER',

	'036':'/DYJetsToLL_M-10To50filter_8TeV-madgraph/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID036_DY10toLLMad-bc58659808d390d23203e20a37543a0e/USER',





# VV

        '074':'/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID074_WZJetsMad_02Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',    
        '075':'/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID075_ZZJetsMad_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',    
        '078':'/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola/amassiro-R52X_S1_V05_S2_V05_S3_V05_ID078_WZTo2L2QMad_13Jun2012-19fe6daf4557731ba3a9d3abfc7c56a5/USER',
        '079':'/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID079_ZZTo2L2QMad_02Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',    
        '082':'/WGstarToLNu2E_TuneZ2star_8TeV-madgraph-tauola/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID082_WGstarToElNuMad_09Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',    
        '083':'/WGstarToLNu2Mu_TuneZ2star_7TeV-madgraph-tauola/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID083_WGstarToMuNuMad_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',    
        '084':'/WGstarToLNu2Tau_TuneZ2star_7TeV-madgraph-tauola/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID084_WGstarToTauNuMad_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',    
        '085':'/WGToLNuG_TuneZ2star_8TeV-madgraph-tauola/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID085_WgammaToLNuG_13Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',    


# DY


	'037c0':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck0-bc58659808d390d23203e20a37543a0e/USER',
        '037c1':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck1-bc58659808d390d23203e20a37543a0e/USER',
        '037c2':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck2-bc58659808d390d23203e20a37543a0e/USER',
        '037c3':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck3-bc58659808d390d23203e20a37543a0e/USER',
        '037c4':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck4-bc58659808d390d23203e20a37543a0e/USER',
        '037c5':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck5-bc58659808d390d23203e20a37543a0e/USER',
        '037c6':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck6-bc58659808d390d23203e20a37543a0e/USER',
        '037c7':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck7-bc58659808d390d23203e20a37543a0e/USER',
        '037c8':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck8-bc58659808d390d23203e20a37543a0e/USER',
        '037c9':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck9-bc58659808d390d23203e20a37543a0e/USER',


# wzttH

        '3110':'/WH_ZH_TTH_HToWW_M-110_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3110_wzttH110ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3115':'/WH_ZH_TTH_HToWW_M-115_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3115_wzttH115ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3118':'wzttH118ToWW',
        '3120':'/WH_ZH_TTH_HToWW_M-120_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3120_wzttH120ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3122':'wzttH122ToWW',
        '3124':'wzttH124ToWW',
        '3125':'/WH_ZH_TTH_HToWW_M-125_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3125_wzttH125ToWW_31May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3126':'wzttH126ToWW',
        '3128':'wzttH128ToWW',
        '3130':'/WH_ZH_TTH_HToWW_M-130_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3130_wzttH130ToWW_02Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3135':'/WH_ZH_TTH_HToWW_M-135_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3135_wzttH135ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3140':'/WH_ZH_TTH_HToWW_M-140_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3140_wzttH140ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3145':'/WH_ZH_TTH_HToWW_M-145_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3145_wzttH145ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3150':'/WH_ZH_TTH_HToWW_M-150_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3150_wzttH150ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3155':'/WH_ZH_TTH_HToWW_M-155_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3155_wzttH155ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3160':'/WH_ZH_TTH_HToWW_M-160_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3160_wzttH160ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3170':'/WH_ZH_TTH_HToWW_M-170_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3170_wzttH170ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3180':'/WH_ZH_TTH_HToWW_M-180_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3180_wzttH180ToWW_31May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3190':'/WH_ZH_TTH_HToWW_M-190_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3190_wzttH190ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3200':'/WH_ZH_TTH_HToWW_M-200_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3200_wzttH200ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3250':'/WH_ZH_TTH_HToWW_M-250_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3250_wzttH250ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3300':'/WH_ZH_TTH_HToWW_M-300_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3300_wzttH300ToWW_31May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3350':'wzttH350ToWW',
        '3400':'wzttH400ToWW',
        '3450':'wzttH450ToWW',
        '3500':'wzttH500ToWW',
        '3550':'wzttH550ToWW',
        '3600':'wzttH600ToWW',   


# VBF

        '2110':'/VBF_HToWWTo2LAndTau2Nu_M-110_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2110_vbfToH110toWWTo2LAndTau2Nu_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2115':'/VBF_HToWWTo2LAndTau2Nu_M-115_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2115_vbfToH115toWWTo2LAndTau2Nu_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2118':'vbfToH118toWWTo2LAndTau2Nu',
        '2120':'/VBF_HToWWTo2LAndTau2Nu_M-120_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2120_vbfToH120toWWTo2LAndTau2Nu_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2122':'vbfToH122toWWTo2LAndTau2Nu',
        '2124':'vbfToH124toWWTo2LAndTau2Nu',
        '2125':'/VBF_HToWWTo2LAndTau2Nu_M-125_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2125_vbfToH125toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2126':'vbfToH126toWWTo2LAndTau2Nu',
        '2128':'vbfToH128toWWTo2LAndTau2Nu',
        '2130':'/VBF_HToWWTo2LAndTau2Nu_M-130_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2130_vbfToH130toWWTo2LAndTau2Nu_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2135':'/VBF_HToWWTo2LAndTau2Nu_M-135_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2135_vbfToH135toWWTo2LAndTau2Nu_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2140':'/VBF_HToWWTo2LAndTau2Nu_M-140_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2140_vbfToH140toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2145':'/VBF_HToWWTo2LAndTau2Nu_M-145_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2145_vbfToH145toWWTo2LAndTau2Nu_02Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2150':'/VBF_HToWWTo2LAndTau2Nu_M-150_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2150_vbfToH150toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2155':'/VBF_HToWWTo2LAndTau2Nu_M-155_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2155_vbfToH155toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2160':'/VBF_HToWWTo2LAndTau2Nu_M-160_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2160_vbfToH160toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2170':'/VBF_HToWWTo2LAndTau2Nu_M-170_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2170_vbfToH170toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2180':'/VBF_HToWWTo2LAndTau2Nu_M-180_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2180_vbfToH180toWWTo2LAndTau2Nu_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2190':'/VBF_HToWWTo2LAndTau2Nu_M-190_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2190_vbfToH190toWWTo2LAndTau2Nu_02Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2200':'/VBF_HToWWTo2LAndTau2Nu_M-200_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2200_vbfToH200toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2250':'/VBF_HToWWTo2LAndTau2Nu_M-250_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2250_vbfToH250toWWTo2LAndTau2Nu_08Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2300':'/VBF_HToWWTo2LAndTau2Nu_M-300_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2300_vbfToH300toWWTo2LAndTau2Nu_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2350':'/VBF_HToWWTo2LAndTau2Nu_M-350_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2350_vbfToH350toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2400':'/VBF_HToWWTo2LAndTau2Nu_M-400_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2400_vbfToH400toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2450':'/VBF_HToWWTo2LAndTau2Nu_M-450_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2450_vbfToH450toWWTo2LAndTau2Nu_02Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2500':'/VBF_HToWWTo2LAndTau2Nu_M-500_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2500_vbfToH500toWWTo2LAndTau2Nu_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2550':'/VBF_HToWWTo2LAndTau2Nu_M-550_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2550_vbfToH550toWWTo2LAndTau2Nu_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2600':'/VBF_HToWWTo2LAndTau2Nu_M-600_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2600_vbfToH600toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2700':'/VBF_HToWWTo2LAndTau2Nu_M-700_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2700_vbfToH700toWWTo2LAndTau2Nu_31May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2800':'/VBF_HToWWTo2LAndTau2Nu_M-800_8TeV-powheg-pythia6/amassiro-R52X_S1_V05_S2_V05_S3_V05_ID2800_vbfToH800toWWTo2LAndTau2Nu_13Jun2012-19fe6daf4557731ba3a9d3abfc7c56a5/USER',
        '2900':'/VBF_HToWWTo2LAndTau2Nu_M-900_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2900_vbfToH900toWWTo2LAndTau2Nu_01Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3000':'/VBF_HToWWTo2LAndTau2Nu_M-1000_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3000_vbfToH1000toWWTo2LAndTau2Nu_07Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',

}
