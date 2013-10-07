
############################
#Make Some Plots
############################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/MakeMELAPlots.C+


############################################
#Create Signal Efficiency
#Need to do it for multiple subchannels
############################################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/yieldAfterMelaCutFit.C+


############################################
#Compute Bkg Yields and Bkg shapes
############################################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/getParams.C+


############################
#Make MELA 2D shape
############################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/create2DHistogram.C+

root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/createPSMelaVsMassHistogram.C+
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/createSpinTwoMinimalMelaVsMassHistogram.C+




############################
#Make Data Cards
############################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/makeSpinCPMeasurementDatacards.C+'(1)'
combineCards.py hzz4l_4mu_7=card_pseudoscalar_2D_m125_7TeV_4mu.txt  hzz4l_4e_7=card_pseudoscalar_2D_m125_7TeV_4e.txt   hzz4l_2e2mu_7=card_pseudoscalar_2D_m125_7TeV_2e2mu.txt  hzz4l_4mu_8=card_pseudoscalar_2D_m125_8TeV_4mu.txt  hzz4l_4e_8=card_pseudoscalar_2D_m125_8TeV_4e.txt     hzz4l_2e2mu_8=card_pseudoscalar_2D_m125_8TeV_2e2mu.txt   >! hzz4lpseudoscalar.text

root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/makeSpinCPMeasurementDatacards.C+'(2)'
combineCards.py hzz4l_4mu_7=card_spinTwoMinimal_2D_m125_7TeV_4mu.txt  hzz4l_4e_7=card_spinTwoMinimal_2D_m125_7TeV_4e.txt   hzz4l_2e2mu_7=card_spinTwoMinimal_2D_m125_7TeV_2e2mu.txt  hzz4l_4mu_8=card_spinTwoMinimal_2D_m125_8TeV_4mu.txt  hzz4l_4e_8=card_spinTwoMinimal_2D_m125_8TeV_4e.txt     hzz4l_2e2mu_8=card_spinTwoMinimal_2D_m125_8TeV_2e2mu.txt   >! hzz4lspinTwoMinimal.text

cp card_pseudoscalar* hzz4lpseudoscalar.text /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v7/
cp card_spinTwoMinimal* hzz4lspinTwoMinimal.text /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/spinTwoMinimal/v7/



############################
#Do Hypothesis Tests
############################

################
#Make Workspace
################
text2workspace.py -m 125 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v7/hzz4lpseudoscalar.text -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o workspace_HCP_pseudoscalar_fixedMu.root
text2workspace.py -m 125 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/spinTwoMinimal/v7/hzz4lspinTwoMinimal.text -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o workspace_HCP_spinTwoMinimal_fixedMu.root

cp workspace_HCP_pseudoscalar_fixedMu.root /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v7/
cp workspace_HCP_spinTwoMinimal_fixedMu.root /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/spinTwoMinimal/v7/


text2workspace.py -m 125 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/card_pseudoscalar_2D_m125_8TeV_combined.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o workspace_HCP_pseudoscalar_fixedMu.root
text2workspace.py -m 125 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/card_pseudoscalar_2D_m125_8TeV_4mu.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o workspace_4mu_HCP_pseudoscalar_fixedMu.root
text2workspace.py -m 125 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/card_pseudoscalar_2D_m125_8TeV_comb.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o workspace_comb_HCP_pseudoscalar_fixedMu.root



#Floating mu
text2workspace.py -m 125 card_pseudoscalar_2D_m125_comb.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs  --PO=muFloating -o workspace_HCP_pseudoscalar_floatMu.root

#2D scan
text2workspace.py -m $MH $card1 -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=muAsPOI -o workspace_HCP_pseudoscalar_twoD.root



################
#Do Toys
################
combine -m 125 -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v6/workspace_HCP_pseudoscalar_fixedMu.root --singlePoint 1 --saveHybridResult --fork 4 -T 4000 -i 1 --clsAcc 0 --fullBToys  


################
#Do data
################
combine -m 125 -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v7/workspace_HCP_pseudoscalar_floatMu.root --singlePoint 1 --saveHybridResult --fork 4 -T 1 -i 1 --clsAcc 0 --fullBToys  
mv higgsCombineTest.HybridNew.mH125.root higgsCombineTest.HybridNew.mH125.pseudoscalar.data.root


combine -m 125 -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/spinTwoMinimal/v7/workspace_HCP_spinTwoMinimal_floatMu.root --singlePoint 1 --saveHybridResult --fork 4 -T 1 -i 1 --clsAcc 0 --fullBToys  
mv higgsCombineTest.HybridNew.mH125.root higgsCombineTest.HybridNew.mH125.spinTwoMinimal.data.root



#Float Mu
combine -M MultiDimFit workspace_HCP_pseudoscalar_floatMu.root --algo=grid --points 100 --fastScan  -m 125 -v 2 -n 1D
combine -M MultiDimFit /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v1/workspace_HCP_pseudoscalar_floatMu.root --algo=grid --points 100 --fastScan  -m 125 -v 2 -n 1D


combine -M MultiDimFit /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v1/workspace_HCP_pseudoscalar_floatMu.root --algo=grid --points 100 --fastScan  -m 125 --setPhysicsModelParameters "x=0.5" -t -1 -v2 -n 1D 



#2D scan
combine -M MultiDimFit /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v1/tworkspace_HCP_pseudoscalar_woD.root --algo=grid --points 10000 --fastScan  -m 125 -v 2


################
#Build qmu from toys
################

foreach f(`seq 1 1 500`)
   root -q -b higgsCombineTest.HybridNew.mH125.${f}.root ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx\(\"qmu.FixedMu.${f}.root\",125,1,\"x\"\)
end

root -q -b higgsCombineTest.HybridNew.mH125.pseudoscalar.data.root ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx\(\"qmu.FixedMu.pseudoscalar.data.root\",125,1,\"x\"\)
cp higgsCombineTest.HybridNew.mH125.pseudoscalar.data.root qmu.FixedMu.pseudoscalar.data.root  /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v7/


root -q -b higgsCombineTest.HybridNew.mH125.spinTwoMinimal.data.root ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx\(\"qmu.FixedMu.spinTwoMinimal.data.root\",125,1,\"x\"\)
cp higgsCombineTest.HybridNew.mH125.spinTwoMinimal.data.root qmu.FixedMu.spinTwoMinimal.data.root  /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/spinTwoMinimal/v7/


################
#Plot Signal Separation
################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/MakePSHiggsHypTestPlots.C+ 

root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/MakeSpinTwoHiggsHypTestPlots.C+ 
