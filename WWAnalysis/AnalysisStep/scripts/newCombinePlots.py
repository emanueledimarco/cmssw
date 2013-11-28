import sys, os, math, re, copy

from WWAnalysis.AnalysisStep.setupSamples_cff import *


# ____             _         __  __       _       
#|  _ \           (_)       |  \/  |     (_)      
#| |_) | ___  __ _ _ _ __   | \  / | __ _ _ _ __  
#|  _ < / _ \/ _` | | '_ \  | |\/| |/ _` | | '_ \ 
#| |_) |  __/ (_| | | | | | | |  | | (_| | | | | |
#|____/ \___|\__, |_|_| |_| |_|  |_|\__,_|_|_| |_|
#             __/ |                               
#            |___/                                


# fileName = sys.argv[1]
# f = TFile(fileName)

c1 = TCanvas("c1","c1")
gStyle.SetOptStat(0)

lumi = 1000.
plots = FWLiteParams.histParams
samplePrefix = "GammaMR"

newChannels = [ chan for chan in channels ]
newChannels.append( "all" )

# masses = [ m for m in allMasses if int(m) == 150 ]
masses = [ "120","160","250"]
for mass in masses:

    f = TFile("data/noteV1L146m{0}.root".format(mass))

#     # make the yield breakdown plots
#     prefix = "eventHists/yields"
#     for chan in newChannels:
#         if chan == "all":
#             hs = getStackedPlot(f,allMassSamples[mass],channels,prefix,"CONVLHT",lumi)
#             addDataPlot(hs,f,channels,prefix,"CONVLHT",-1)
#         else:
#             hs = getStackedPlot(f,allMassSamples[mass],[chan],prefix,"CONVLHT",lumi)
#             addDataPlot(hs,f,[chan],prefix,"CONVLHT",-1)
# 
#         c1.SetLogy(False)
#         hist = hs.Draw(1)
#         c1.Print("plots/"+chan+"/"+plotName+"/{0:0>2.0f}".format(cut)+"_"+fullSet[cut].simple.value()+".lin.png")
#         c1.Print("plots/"+chan+"/"+plotName+"/{0:0>2.0f}".format(cut)+"_"+fullSet[cut].simple.value()+".lin.eps")
#         
#         c1.SetLogy(True)
#         hist = hs.Draw(1)
#         c1.Print("plots/"+chan+"/"+plotName+"/{0:0>2.0f}".format(cut)+"_"+fullSet[cut].simple.value()+".log.png")
#         c1.Print("plots/"+chan+"/"+plotName+"/{0:0>2.0f}".format(cut)+"_"+fullSet[cut].simple.value()+".log.eps")


    # Make the kine plots
    prefix = "eventHists/bycut"
    fullSet = cloneVPSet(defaultWW)
    addMassDependentCuts(fullSet,getattr(sys.modules[__name__],"h{0}{1}".format(samplePrefix,mass)))
    # injectTightIsolationCut(fullSet)
    injectIPCut(fullSet)

    for plotName in plots.parameterNames_():
        plotInfo = getattr(plots,plotName)

#         for cut in range(0,len(fullSet)):
        for cut in [12]:
            postfix = plotName+"/"+"{0:0>2.0f}".format(cut)

            for chan in newChannels:
                if chan == "all":
                    hs = getStackedPlot(f,allMassSamples[mass],channels,prefix,"CONVLHT/"+postfix,lumi)
#                     addDataPlot(hs,f,channels,prefix,"CONVLHT/"+postfix,-1)
                else:
                    continue
                    hs = getStackedPlot(f,allMassSamples[mass],[chan],prefix,"CONVLHT/"+postfix,lumi)
#                     addDataPlot(hs,f,[chan],prefix,"CONVLHT/"+postfix,-1)
                hs.addLabel("After " + fullSet[cut].label.value())
                hs.setMass(int(mass))

#                 rebin = 1 if plotInfo.nbins.value() < 100 else 5
                rebinTo = 0 if plotInfo.nbins.value() < 25 else 20
                c1.SetLogy(False)
#                 hist = hs.Draw(rebin)
                hist = hs.DrawAndRebinTo(rebinTo)
                c1.Print("plots/"+chan+"/"+plotName+"/{0:0>2.0f}".format(cut)+"_"+fullSet[cut].simple.value()+"_{0}_lin.png".format(mass))
                c1.Print("plots/"+chan+"/"+plotName+"/{0:0>2.0f}".format(cut)+"_"+fullSet[cut].simple.value()+"_{0}_lin.eps".format(mass))
                
                c1.SetLogy(True)
                hist = hs.Draw(1)
                c1.Print("plots/"+chan+"/"+plotName+"/{0:0>2.0f}".format(cut)+"_"+fullSet[cut].simple.value()+"_{0}_log.png".format(mass))
                c1.Print("plots/"+chan+"/"+plotName+"/{0:0>2.0f}".format(cut)+"_"+fullSet[cut].simple.value()+"_{0}_log.eps".format(mass))



