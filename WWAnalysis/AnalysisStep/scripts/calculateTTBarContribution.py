
from ROOT import gROOT, TCanvas, TH1F, TFile, TDirectoryFile, TLegend, gInterpreter
from ROOT import kCyan, kGreen, kMagenta, kBlue, kBlack, kRed, gStyle
from ROOT import TEfficiency, Math
from WWAnalysis.AnalysisStep.scaleFactors_cff import ttbarSamples, dataSamples

import sys, math, os

gInterpreter.ExecuteMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoStyle.C" )

def myFormat(v,e):
    if e>0:
        le = int(math.floor(math.log(e,10)))+1
        right = -1*min(0,le-2)
    else:
        right = 1
    tString = '{0:.'+'{0}'.format(right)+'f},{1:.'+'{0}'.format(right)+'f}'
    return tString.format(v,e)


f = TFile(sys.argv[1])
yieldDir = f.Get('eventHists/yields')

channels = ['wwmumu','wwmuel','wwelmu','wwelel']
types =                 ['TopBMC','TopMuMC', 'TopTagMC']
labels = dict(zip(types,['b-tag', 'soft-#mu','combined']))
colors =                [ kBlue,  kRed      , kBlack]

# mcHists = "023.TTJetsMad"
# hist = "101160.ggToH160toWWto2L2Nu"

# c1 = TCanvas()
# gStyle.SetOptStat(0)

# First make the SMURF plot
# numHists = {}
# denHists = {}
# divHists = {}
# drawString = ""
# l = TLegend(0.65,0.55,0.9,0.8)
# l.SetBorderSize(     0);
# l.SetFillColor (     0);
# l.SetTextAlign (    12);
# l.SetTextFont  (    42);
# l.SetTextSize  (  0.03);
# l.SetHeader("t#bar{t} MC")
# l.SetFillColor(0)
# 
# for type in types:
#     for chan in channels:
#         for sample in ttbarSamples:
#             thisNum = "eventHists/bycut/"+chan+type+"/nJets/08/"+sample+"."+ttbarSamples[sample][0]
#             thisDen = "eventHists/bycut/"+chan+type+"/nJets/07/"+sample+"."+ttbarSamples[sample][0]
# 
#             h = f.Get(thisNum)
#             if type in numHists:
#                 numHists[type].Add(h, ttbarSamples[sample][1])
#             else: 
#                 numHists[type] = h.Clone()
#                 numHists[type].Scale(ttbarSamples[sample][1])
# 
#             h = f.Get(thisDen)
#             if type in denHists:
#                 denHists[type].Add(h, ttbarSamples[sample][1])
#             else: 
#                 denHists[type] = h.Clone()
#                 denHists[type].Scale(ttbarSamples[sample][1])
#     
#     divHists[type] = numHists[type].Clone( "div"+type )
#     divHists[type].Divide(numHists[type],denHists[type])
#     divHists[type].SetLineColor(colors[-1])
#     divHists[type].SetLineWidth(2)
#     divHists[type].Draw(drawString)
#     divHists[type].GetXaxis().SetRangeUser(0,5)
#     divHists[type].SetMinimum(0.0)
#     divHists[type].SetMaximum(1.2)
#     divHists[type].GetYaxis().SetTitle("Efficiency")
#     divHists[type].GetXaxis().SetNdivisions(106)
#     divHists[type].SetTitle("CMS preliminary")
#     l.AddEntry(divHists[type],labels[type],"l")
#     drawString = "SAME"
#     colors.pop()
# 
# l.Draw()
# c1.Print("topTagging.smurfFig20.png")
# c1.Print("topTagging.smurfFig20.eps")

sigma = 0.68540158589942957
alpha = (1.0 - sigma)/2.

def getEffandError(num,den):
    eff = 0 if den == 0 else num/den
    hi = TEfficiency.ClopperPearson(int(den),int(num),sigma,True)
    lo  = TEfficiency.ClopperPearson(int(den),int(num),sigma,False)
    return (eff, (eff-lo)/eff, (hi-eff)/eff)

def getEffFromHist(type,samples,addWeight):

    if not addWeight:
        num = den = 0
        for sample in samples:
            for chan in channels:
                h = yieldDir.Get(chan+type).Get(sample+"."+samples[sample][0])
                num += h.GetBinContent(h.GetNbinsX())
                den += h.GetBinContent(h.GetNbinsX()-1)
        print num,den,getEffandError(num,den)
        return getEffandError(num,den)

    effs = {}
    for sample in samples:
        num = den = 0
        for chan in channels:
            h = yieldDir.Get(chan+type).Get(sample+"."+samples[sample][0])
            num += h.GetBinContent(h.GetNbinsX())
            den += h.GetBinContent(h.GetNbinsX()-1)
        effs[sample] = getEffandError(num,den) 

    weight = eff = lo = hi = 0
    for sample in samples:
        weight += samples[sample][1]
        eff += effs[sample][0] * samples[sample][1] 
        lo  += effs[sample][1] * effs[sample][1] * samples[sample][1] * samples[sample][1] 
        hi  += effs[sample][2] * effs[sample][2] * samples[sample][1] * samples[sample][1] 
    eff /= weight
    lo = math.sqrt(lo) / weight
    hi = math.sqrt(hi) / weight
    return (eff,lo,hi)
    

# eff1bMC     = getEffFromHist('TopB'  , ttbarSamples, True) 
eff1bData   = getEffFromHist('TopB'  , dataSamples, False)
# eff2bMC     = getEffandError(numHists['TopBMC'  ].GetBinContent(1),denHists['TopBMC'  ].GetBinContent(1))
# effSoft     = getEffandError(numHists['TopMuMC' ].GetBinContent(1),denHists['TopMuMC' ].GetBinContent(1))
# eff2bFrom1b = 1 - (1 - eff1bMC[0]  ) * (1 - eff1bMC[0]   )
# eff2bData   = 1 - (1 - eff1bData[0]) * (1 - eff1bData[0] )
# 
# print "\\begin{tabular}{|l|c|} \hline"
# print "  Actual MC b-tag ($\epsilon_{{1b}}^{{mc}}$):      & ${0:.2%}_{{ -{1:.2%} }}^{{ +{2:.2%} }}$ \\\\ [1.0mm]".format(eff2bMC[0],eff2bMC[0]*eff2bMC[1],eff2bMC[0]*eff2bMC[2])
# print "  1-jet MC b-tag ($\epsilon_{{2b}}^{{mc}}$):       & ${0:.2%}_{{ -{1:.2%} }}^{{ +{2:.2%} }}$ \\\\ [1.0mm]".format(eff1bMC[0],eff1bMC[0]*eff1bMC[1],eff1bMC[0]*eff1bMC[2])
# print "  1-jet data b-tag ($\epsilon_{{2b}}^{{data}}$):   & ${0:.2%}_{{ -{1:.2%} }}^{{ +{2:.2%} }}$ \\\\ [1.0mm]".format(eff1bData[0],eff1bData[0]*eff1bData[1],eff1bData[0]*eff1bData[2])
# print "  Actual MC soft-$\mu$ ($\epsilon_{{2b}}^{{mc}}$): & ${0:.2%}_{{ -{1:.2%} }}^{{ +{2:.2%} }}$ \\\\ [1.0mm]".format(effSoft[0],effSoft[0]*effSoft[1],effSoft[0]*effSoft[2]),
# print "\hline"
# 
# print "  b-tag closure test:                              & ${0:.2%}_{{ -{1:.2%} }}^{{ +{2:.2%} }}$ \\\\ [1.0mm]".format(eff2bFrom1b,0,0)
# print "  Efficiency from data ($\epsilon_{{2b}}^{{data}}$):   & ${0:.2%}_{{ -{1:.2%} }}^{{ +{2:.2%} }}$ \\\\ [1.0mm]".format(eff2bData,0,0),
# print "\hline"
# print "\\end{tabular}"
# print

# Also want to eventually make the purity plot

