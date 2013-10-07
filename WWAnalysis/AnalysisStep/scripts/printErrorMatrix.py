from ROOT import gROOT, TCanvas, TH1F, TFile, TDirectoryFile
from ROOT import kCyan, kGreen, kMagenta, kBlue

import sys
import math 

def myFormat(v,e):
    if e>0:
        le = int(math.floor(math.log(e,10)))+1
        right = -1*min(0,le-2)
    else:
        right = 1
    tString = '{0:.'+'{0}'.format(right)+'f},{1:.'+'{0}'.format(right)+'f}'
    return tString.format(v,e)


if len(sys.argv) < 2:
    print "I need a file!"
    sys.exit()

lumi=35.49
if(len(sys.argv)>=3):
    lumi=float(sys.argv[2])

from WWAnalysis.AnalysisStep.data.v03.higgsSamples_cff import higgsSamples
from WWAnalysis.AnalysisStep.data.v03.failureSamples_cff import failureSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff import wwSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff import wwGroupings
from WWAnalysis.AnalysisStep.data.v03.sep172010BSamples_cff import sep172010BSamples

f = TFile(sys.argv[1])
scaledDir = f.Get('scaled')
rawDir = f.Get('raw')
sampleKeys = [key.GetName() for key in scaledDir.GetListOfKeys()]

nbins = 1
channels = ['wwelel','wwelmu','wwmumu','all']

backgrounds = [ wwSamples[k][0] for k in wwSamples ]
signals = [ higgsSamples[k][0] for k in higgsSamples ] 
failures = [ failureSamples[k][0] for k in failureSamples ] 
dataSamples = [ sep172010BSamples[k][0] for k in sep172010BSamples ]
allsamples = backgrounds + signals + dataSamples + failures

enhancements =      dict(zip(signals, [ higgsSamples[k][4] for k in higgsSamples ]))
enhancements.update(dict(zip(failures,[ failureSamples[k][4] for k in failureSamples ])))

scales =       dict(zip(signals,     [ higgsSamples[k][1]       for k in higgsSamples ]))
scales.update( dict(zip(failures,     [ failureSamples[k][1]     for k in failureSamples     ])))
scales.update( dict(zip(dataSamples, [ sep172010BSamples[k][1]  for k in sep172010BSamples  ])))
scales.update( dict(zip(backgrounds, [ wwSamples[k][1]          for k in wwSamples          ])))

allGroupings = {}
allGroupings["hww"] = [val for val in signals if val in [ k.replace("all_",'',1) for k in sampleKeys if 'all_' in k]]
allGroupings.update( wwGroupings )
allGroupings["data"] = [val for val in dataSamples if val in [ k.replace("all_",'',1) for k in sampleKeys if 'all_' in k]]

bkgYields = {}
sigYields = {}
enYields = {}
dataYields = {}
bkgErr2 = {}
sigErr2 = {}
enErr2 = {}
dataErr2 = {}

countsScaled = {}
errorsSquaredScaled = {}

for channel in channels:
    countsScaled[channel] = {}
    errorsSquaredScaled[channel] = {}
    for group in allGroupings:
        countsScaled[channel][group] = 0.0
        errorsSquaredScaled[channel][group] = 0.0

    for k in sampleKeys:
        if k.split('_')[0] in channel:
            countsScaled[channel][k] = {}
            errorsSquaredScaled[channel][k] = {}
            hScaled = scaledDir.Get(k)
            hRaw = rawDir.Get(k)
            sample = k.replace(channel+'_','',1)
            for group in allGroupings:
                if sample in allGroupings[group]:
                    countsScaled[channel][group] += hScaled.GetBinContent(1)*scales[sample]*lumi
                    errorsSquaredScaled[channel][group] += pow( hScaled.GetBinError(1)*scales[sample]*lumi, 2)

                

print 'Labels'.center(30),
for group in ['hww','qqWW','ggWW','diBoson','top','zJets','wJets']:
    if countsScaled["all"][group] != 0:
        print group.center(10),
print

print '{0} Yield'.format(allGroupings['hww'][0][-3:]).center(30),
for group in ['hww','qqWW','ggWW','diBoson','top','zJets','wJets']:
    if countsScaled["all"][group] != 0:
        print myFormat(countsScaled["all"][group],math.sqrt(errorsSquaredScaled["all"][group])/countsScaled["all"][group]).split(',')[0].center(10),
print

print '{0} Error'.format(allGroupings['hww'][0][-3:]).center(30),
for group in ['hww','qqWW','ggWW','diBoson','top','zJets','wJets']:
    if countsScaled["all"][group] != 0:
        print myFormat(countsScaled["all"][group],math.sqrt(errorsSquaredScaled["all"][group])/countsScaled["all"][group]).split(',')[1].center(10),
print

