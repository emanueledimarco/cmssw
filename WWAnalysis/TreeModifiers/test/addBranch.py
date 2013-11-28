import FWCore.ParameterSet.Config as cms
import math

process = cms.Process("BranchAdded")

process.General = cms.PSet(   
    inputFiles    =     cms.vstring('file:hzzTree_id113.root'),
    branchAdders    =     cms.vstring(#'BaseWeightBranchAdder',
                                      'PileupBranchAdder',
                                      #'LeptSfBranchAdder',
                                      #'LeptPRBranchAdder',
                                      ),


)


process.BaseWeightBranchAdder = cms.PSet(
    treesToProcess = cms.vstring(), # if empty process all the trees in the file
    newBranchName = cms.string("baseW"),
    datasetConf = cms.string("8TeV"),
)


process.PileupBranchAdder = cms.PSet(
    treesToProcess = cms.vstring('zz4lTree','zxTree','anyZxTree',
                                 'zmmtree','zeetree',
                                 'zllmtree','zlletree'),
    newBranchName = cms.string("puW"),
    pileupConf = cms.string("2012"),
)

process.LeptSfBranchAdder = cms.PSet(
    treesToProcess = cms.vstring('zz4lTree','zxTree','anyZxTree'),
    newBranchName = cms.string("lepW"),
    datasetConf = cms.string("2012"),
)

process.LeptPRBranchAdder = cms.PSet(
    treesToProcess = cms.vstring('zxTree'),
    newBranchName = cms.string("lepPRW"),
    datasetConf = cms.string("2012"),
)

### HERE YOU HAVE TO ADD THE CONFIGURATION FOR ANY NEW BRANCHADDER
#process.myBranchAdder = cms.PSet(
#    treesToProcess = cms.vstring('zz4lTree'),
#    newBranchName = cms.string("puW"),
#)



#############################
#CUSTOMIZATION
#############################

process.BaseWeightBranchAdder.datasetConf = cms.string("8TeV")

process.PileupBranchAdder.pileupConf = cms.string("puV10toIchep53X")
process.PileupBranchAdder.newBranchName = cms.string("puICHEP")
#process.PileupBranchAdder.pileupConf = cms.string("puV10toHcp53X")
#process.PileupBranchAdder.newBranchName = cms.string("puHCP")

process.LeptSfBranchAdder.datasetConf = cms.string("2012")
process.LeptPRBranchAdder.datasetConf = cms.string("2012")

process.General.inputFiles = cms.vstring(
"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id3275.root",
"file:/dummy"
)

#ls -lhtr /data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02/MC/8TeV/|awk '{print "\"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02/MC/8TeV/"$9"\","}'
#process.General.inputFiles = cms.vstring(
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id110.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id111.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id114.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id115.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id116.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id117.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id118.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id119.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id120.root",
###
##ggH signals
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1115.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1116.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1117.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1118.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1119.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1120.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1121.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1122.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1124.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1125.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1126.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1127.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1129.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1130.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1135.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1140.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1145.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1150.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1160.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1170.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1180.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1190.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1220.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1250.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1275.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1300.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1325.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1350.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1375.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1400.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1425.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1450.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1475.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1500.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1525.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1550.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1600.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1650.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1700.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1750.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1800.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1850.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1900.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id1950.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id11000.root",
###
###VBF signals
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2115.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2116.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2117.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2118.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2119.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2120.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2121.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2122.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2123.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2124.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2125.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2126.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2127.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2128.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2129.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2130.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2135.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2140.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2145.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2150.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2160.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2170.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2180.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2190.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2200.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2220.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2250.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2275.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2325.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2350.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2375.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2400.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2425.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2450.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2500.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2525.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2550.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2575.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2600.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2650.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2700.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2750.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2800.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2850.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2900.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id2950.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id21000.root",
###
###associated-production signals
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3115.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3116.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3117.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3118.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3119.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3120.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3121.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3122.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3123.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3124.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3125.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3126.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3127.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3128.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3129.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3130.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3135.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3140.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3145.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3150.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3160.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3170.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3180.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3190.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3200.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3220.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3250.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id3275.root"
#)
#
# pseudo-scalar higgs
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id81000.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id8110.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id8120.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id8125.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id8130.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id8140.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id8250.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id8500.root",
#
# graviton
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id9110.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id9120.root",
#"file:/data/hzz4l/step3/HZZ4L_53X_S1_V11_S2_V02_S3_V00/MC/8TeV/noRegrYesCalibYesMu/hzzTree_id9125.root",





#ls -lhtr /data/hzz4l/step3/2011/MC/HZZ4L_S3_V01/* |awk '{print "\"file:/data/hzz4l/step3/2011/MC/HZZ4L_S3_V01/"$9"\","}'

#process.BaseWeightBranchAdder.datasetConf = cms.string("7TeV")
#process.PileupBranchAdder.pileupConf = cms.string("2011")
#process.LeptSfBranchAdder.datasetConf = cms.string("2011")
#process.LeptPRBranchAdder.datasetConf = cms.string("2011")
#process.General.inputFiles = cms.vstring(
#signal files
#"file:/data/hzz4l/step3/2011/MC/HZZ4L_S3_V02/hzzTree_id200.root",
