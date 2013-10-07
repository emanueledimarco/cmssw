import sys, os, math, re, copy

from ROOT import gROOT, gStyle, TCanvas, TH1F, TFile, TDirectoryFile, TStyle, TLegend, THStack
from ROOT import kCyan, kGreen, kMagenta, kBlue, kWhite, kOrange, kViolet, kBlack
from ROOT import TObject, gInterpreter, Long
from WWAnalysis.AnalysisStep.scaleFactors_cff import *
from WWAnalysis.AnalysisStep.yieldProducer_cfi import FWLiteParams
from WWAnalysis.AnalysisStep.cutPSets_cfi import *
from WWAnalysis.AnalysisStep.cutMassDependPSets_cfi import *


gROOT.LoadMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoPlot.C+" )
gInterpreter.ExecuteMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoStyle.C" )
from ROOT import LatinoPlot


# # loopables and other default that can probably be configured later
# allMasses = [ str(x) for x in [120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600] ]
# 
# # the suffix of whatever the name of the SkimEvent branch was, i.e. wwelelIPLHT
# suffixes = ['CONVLHT']
# 
# # four channels plus their corresponding LaTeX labels
# channels =                        ['wwmumu'  ,'wwmuel','wwelmu','wwelel'  ]
# channelLabels = dict(zip(channels,['$\mu\mu$','$\mu$e','e$\mu$','ee']))
# 
# # Sample setup
# # sampleLabels =                       [      'WW',        'WZ/ZZ',  'Z+Jets',        'Top',    'W+Jets', ]
# # myEnum       = dict(zip(sampleLabels,[         1,              4,         2,            3,           5, ]))
# # tempSamples  = dict(zip(sampleLabels,[ wwSamples, diBosonSamples, dySamples, ttbarSamples, wJetSamples, ]))
# sampleLabels =                       [       'HWW',      'WW',        'WZ/ZZ',  'Z+Jets',        'Top',    'W+Jets', ]
# myEnum       = dict(zip(sampleLabels,[           0,         1,              4,         2,            3,           5, ]))
# tempSamples  = dict(zip(sampleLabels,[ h160Samples, wwSamples, diBosonSamples, dySamples,   topSamples, wJetSamples, ]))
# 
# allMassSamples = {}
# for mass in allMasses:
#     allMassSamples[mass] = copy.deepcopy(tempSamples)
#     allMassSamples[mass]['HWW'] = getattr(sys.modules[__name__],"h{0}Samples".format(mass))
# 
# 
# # Sample setup
# # cardChannels = ['elel'  ,'mumu'  ,'elmu'  ,'muel'  ]
# cardChannels = ['mumu'  ,'muel'  ,'elmu' , 'elel' ]
# cardLabels =                           [       'HWW',     'WW',     'ggWW',           'VV',        'Top',    'Zjet',      'Wjet',]  
# tempCardSamples  = dict(zip(cardLabels,[ h160Samples, wwSamples, ggWWSample, diBosonSamples,   topSamples, dySamples, wJetSamples,]))
# isDataEstimate   = dict(zip(cardLabels,[ False      , True    , False     , False         ,   True      , True     , True       ,]))
# # cardLabels =                           [       'HWW',     'WW',     'ggWW',           'VV',        'Top',    'Zjet',      'Wjet',      'Wgam', ]
# # tempCardSamples  = dict(zip(cardLabels,[ h160Samples, wwSamples, ggWWSample, diBosonSamples,   topSamples, dySamples, wJetSamples, wGamSamples, ]))
# 
# cardMassSamples = {}
# for mass in allMasses:
#     cardMassSamples[mass] = copy.deepcopy(tempCardSamples)
#     cardMassSamples[mass]['HWW'] = getattr(sys.modules[__name__],"h{0}Samples".format(mass))
# 
# # Sample setup
# 
# latexLabels =                            [       'HWW',      'WW',           'WZ',           'ZZ',  'Z+Jets',         'tt',   'single t',    'W+Jets', ]
# tempLatexSamples  = dict(zip(latexLabels,[ h160Samples, wwSamples,      wzSamples,      zzSamples, dySamples, ttbarSamples,    tWSamples, wJetSamples, ]))
# 
# latexMassSamples = {}
# for mass in allMasses:
#     latexMassSamples[mass] = copy.deepcopy(tempLatexSamples)
#     latexMassSamples[mass]['HWW'] = getattr(sys.modules[__name__],"h{0}Samples".format(mass))
# 
# bkgNoteLabels =                             [   'Z+jets',  '$\\mathrm{t}\\bar{\\mathrm{t}}$',  'single t',    'W+jets',        'WZ/ZZ',     'ggWW',   'qqWW', ]
# tempBkgNoteSamples = dict(zip(bkgNoteLabels,[  dySamples, ttbarSamples,   tWSamples, wJetSamples, diBosonSamples, ggWWSample, wwSamples, ]))
# 
# bkgNoteMassSamples = {}
# for mass in allMasses:
#     bkgNoteMassSamples[mass] = copy.deepcopy(tempBkgNoteSamples)
#     bkgNoteMassSamples[mass]['HWW'] = getattr(sys.modules[__name__],"h{0}Samples".format(mass))
# 
# def getCombinedHistogram(f,sample,channels,prefix,suffix,lumi=1000):
#     returnPlot = None
#     for channel in channels:
#         for key in sample:
#             plotName = prefix+"/"+channel+suffix+"/"+key+"_"+sample[key][0]
#             thisPlot = f.Get(plotName)
# 
#             if thisPlot == None: 
#                 print "WARNING: {0} wasn't found in {1}".format(plotName,f.GetPath())
#                 continue
# 
#             if returnPlot != None:
#                 if lumi < 0: returnPlot.Add(thisPlot)
#                 else:        returnPlot.Add(thisPlot,lumi/1000*sample[key][1])
#             else:
#                 returnPlot = thisPlot.Clone(key+"stack")
#                 if lumi > 0: returnPlot.Scale(lumi/1000*sample[key][1])
#     return returnPlot
#     
# def getStackedYieldPlot(f,samples,channels,lumi=1000):
#     addablePlots = {}
#     for label in samples:
#         addablePlots[label] = getCombinedHistogram(f,samples[label],channels,prefix,suffix,lumi)
#     
#     hs = LatinoPlot()
#     hs.setLumi(lumi)
#     re1='(.*?)' # should be the label
#     re2='(\\[)' # left bracket
#     re3='(.*?)' # units
#     re4='(\\])' # right bracket
#     rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
#     for label in samples:
#         if addablePlots[label] == None: continue
#         hs.setMCHist(Long(myEnum[label]),addablePlots[label])
#         # I hate myself for this ...
#         m = rg.search(addablePlots[label].GetXaxis().GetTitle())
#         if m:
#             hs.setLabel(m.group(1))
#             hs.setUnits(m.group(3))
#         elif addablePlots[label].GetXaxis().GetBinLabel(1) == "":
#             hs.setLabel(addablePlots[label].GetXaxis().GetTitle())
#             hs.setUnits("")
#         else:
#             hs.setBreakdown()
#     return hs
#     
# def getStackedPlot(f,samples,channels,prefix,suffix,lumi=1000):
#     addablePlots = {}
#     for label in samples:
#         addablePlots[label] = getCombinedHistogram(f,samples[label],channels,prefix,suffix,lumi)
#     
#     hs = LatinoPlot()
#     hs.setLumi(lumi)
#     re1='(.*?)' # should be the label
#     re2='(\\[)' # left bracket
#     re3='(.*?)' # units
#     re4='(\\])' # right bracket
#     rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
#     for label in samples:
#         if addablePlots[label] == None: continue
#         hs.setMCHist(Long(myEnum[label]),addablePlots[label])
#         # I hate myself for this ...
#         m = rg.search(addablePlots[label].GetXaxis().GetTitle())
#         if m:
#             hs.setLabel(m.group(1))
#             hs.setUnits(m.group(3))
#         elif addablePlots[label].GetXaxis().GetBinLabel(1) == "":
#             hs.setLabel(addablePlots[label].GetXaxis().GetTitle())
#             hs.setUnits("")
#         else:
#             hs.setBreakdown()
#     return hs
#     
# def addDataPlot(hs,f,channels,prefix,suffix,lumi):
# #     combinedData = None
# #     for key in dataSamples:
# #         thisPlot = f.Get(plotPath+key+"."+dataSamples[key][0])
#     combinedData = getCombinedHistogram(f,dataSamples,channels,prefix,suffix,lumi)
#     hs.setDataHist(combinedData)
#     
# 
# def getYieldAndError(f,sample,channels,cut=-1,suffix=suffixes[0],lumi=1000):
#     prefix = "eventHists/yields"
#     combinedHist = getCombinedHistogram(f,sample,channels,prefix,suffix,lumi)
#     bin = cut if cut != -1 else combinedHist.GetNbinsX()
#     if combinedHist == None:
#         return (0,0)
#     else:
#         return (combinedHist.GetBinContent(bin),combinedHist.GetBinError(bin))
# 
# def getYield(f,sample,channels,cut=-1,suffix=suffixes[0],lumi=1000):
#     return getYieldAndError(f,sample,channels,cut,suffix,lumi)[0]
# 
# def getError(f,sample,channels,cut=-1,suffix=suffixes[0],lumi=1000):
#     return getYieldAndError(f,sample,channels,cut,suffix,lumi)[1]
# 
# def convertToRealLatex(s):
# 
#     s = s.replace("#slash{E}_{T}","\\met")
#     s = s.replace("<","$<$")
#     s = s.replace(">","$>$")
# 
#     txt='#adswf'
#     
#     re1='(#)'   # Any Single Character 1
#     re2='((?:[a-z][a-z]+))' # Word 1
#     
#     rg = re.compile(re1+re2,re.IGNORECASE|re.DOTALL)
#     m = rg.search(s)
#     while m:
#         c1=m.group(1)
#         word1=m.group(2)
#         s = s.replace(c1+word1,"$\\"+word1+"$")
#         m = rg.search(s)
# 
# #     re1='(((?<!\$)?:[a-z]+' # Word 1
#     re1='((?:[a-z_\^\{\}\*\/]+' # Word 1
#     re2='[_\^]'   # Any Single Character 1
#     re3='\\{+.*?\\}+))' # Word 1
#     
#     rg = re.compile(re1+re2+re3,re.IGNORECASE|re.DOTALL)
#     m = rg.search(s)
#     if m:
#         word=m.group(1)
#         s = s.replace(word,"$"+word+"$")
#     
# 
#     return s
# 
# 
