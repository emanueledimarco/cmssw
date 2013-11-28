import sys, os, math, re, copy

from WWAnalysis.AnalysisStep.setupSamples_cff import *

def getNuissances(fn):
    f = open(fn,'r')
    lines = f.readlines()

    channels = lines[0].split()
    samples  = lines[1].split()
    if len(channels) != len(samples):
        print "ERROR, bad freakin nuissance file dude.  Just copy it from the damn spreadsheet, jeeez."
        sys.exit(1)

    # build the 3D matrix: channel x sample x nuissances
    theNuis = {}
    for chan in channels:
        theNuis[chan] = {}

    labels = []
    types = []
    norms = []
    for line in lines[2:]:
        if len(line.split()[3:]) != len(channels):
            print "ERROR, bad freakin nuissance file dude.  Just copy it from the damn spreadsheet, jeeez."
            sys.exit(1)

        columns = line.split()
        labels.append(columns[0])
        types.append(columns[1])
        norms.append(columns[2])
        nuises = columns[3:]

        for i in range(len(channels)):
            if samples[i] in theNuis[channels[i]]:
                theNuis[channels[i]][samples[i]].append( nuises[i] )
            else:
                theNuis[channels[i]][samples[i]] = [ nuises[i] ]

    return (labels,types,norms,theNuis)

def printDashString(outfile,a,b,value):
    for i in range(len(cardChannels)):
        for j in range(len(cardLabels)):
            if i == a and j == b: print >> outfile, value,
            else:                 print >> outfile, "{0:^8}".format('-'),
    print >> outfile
# ____             _         __  __       _       
#|  _ \           (_)       |  \/  |     (_)      
#| |_) | ___  __ _ _ _ __   | \  / | __ _ _ _ __  
#|  _ < / _ \/ _` | | '_ \  | |\/| |/ _` | | '_ \ 
#| |_) |  __/ (_| | | | | | | |  | | (_| | | | | |
#|____/ \___|\__, |_|_| |_| |_|  |_|\__,_|_|_| |_|
#             __/ |                               
#            |___/                                

import commands
def getDataEstimateAndError(mass,chan,samp,jet):
    est = commands.getoutput('grep {0} {1}.{2}.txt | c {3}'.format(mass,samp,jet,cardChannels.index(chan)*2+2))
    err = commands.getoutput('grep {0} {1}.{2}.txt | c {3}'.format(mass,samp,jet,cardChannels.index(chan)*2+3))
    return (0 if est == "" else float(est),0 if err == "" else math.fabs(float(err)))

def getDataEstimate(mass,chan,samp,jet):
    return getDataEstimateAndError(mass,chan,samp,jet)[0]

def getDataError(mass,chan,samp,jet):
    return getDataEstimateAndError(mass,chan,samp,jet)[1]

# assumes file is like our spreadsheet
(nuisLabels,nuisTypes,norms,nuisValues) = getNuissances('nuissances.txt')


for mass in allMasses:
#     f = TFile('workingDir/hists/data/note{0}.root'.format(mass))
    f = TFile('../workingDir/hists/data/noteV1L146m{0}.root'.format(mass))
#     outFile = open('cards/wVGamma/ww.{0}.txt'.format(mass),'w')
    outFile = open('cards/ww.{0}.txt'.format(mass),'w')
    
    # header
    print >> outFile, "imax {0}".format(len(nuisValues))
    print >> outFile, "jmax *"
    print >> outFile, "kmax *"
    print >> outFile, "------------"
    
    # observation
    print >> outFile, "{0:<24}".format("bin"),
    for chan in cardChannels:
        print >> outFile, "{0:^8}".format(chan),
    print >> outFile
    
    print >> outFile, "{0:<24}".format("observation"),
    for chan in cardChannels:
        sum = 0
        # not the higgs
        for samp in cardLabels[1:]:
            sum += getYield(f,cardMassSamples[mass][samp],['ww'+chan],-1,"CONVLHT")
        print >> outFile, "{0:^8}".format(int(math.floor(sum))),
    print >> outFile
    print >> outFile, "------------"
    
    # yields
    print >> outFile, "{0:<24}".format("bin"),
    for chan in cardChannels:
        for samp in cardLabels:
            print >> outFile, "{0:^8}".format(chan),
    print >> outFile
    
    print >> outFile, "{0:<24}".format("process"),
    for chan in cardChannels:
        for samp in cardLabels:
            print >> outFile, "{0:^8}".format(samp),
    print >> outFile
    
    print >> outFile, "{0:<24}".format("process"),
    for chan in cardChannels:
        for i in range(len(cardLabels)):
            print >> outFile, "{0:^8}".format(i),
    print >> outFile
    
    print >> outFile, "{0:<24}".format("rate"),
    for chan in cardChannels:
        for samp in cardLabels:
            if isDataEstimate[samp]:
                print >> outFile, "{0:^8.3f}".format( getDataEstimate(mass,chan,samp,0) ),
            else:
                print >> outFile, "{0:^8.3f}".format( getYield(f,cardMassSamples[mass][samp],['ww'+chan],-1,"CONVLHT") ),
    print >> outFile
    print >> outFile, "------------"
    
    # fill me with MC (or data ??) stats
    for i in range(len(cardChannels)):
        for j in range(len(cardLabels)):
            print >> outFile, "{0:<14}{1:<10}".format("stat"+cardLabels[j]+cardChannels[i],"lnN"),
            yieldAndError = ()
            if isDataEstimate[cardLabels[j]]:
                yieldAndError = getDataEstimateAndError(mass,cardChannels[i],cardLabels[j],0)
            else:
                yieldAndError = getYieldAndError(f,cardMassSamples[mass][cardLabels[j]],['ww'+cardChannels[i]],-1,"CONVLHT")
            printDashString(outFile,i,j,"{0:^8.3f}".format( 1. + yieldAndError[1]/yieldAndError[0] if yieldAndError[0] != 0 else 1.0))
    
    for i in range(len(nuisLabels)):
        if nuisTypes[i] == "lnN":
            print >> outFile, "{0:<14}{1:<10}".format(nuisLabels[i],nuisTypes[i]),
            for chan in cardChannels:
                for samp in cardLabels:
                    print >> outFile, "{0:^8}".format(nuisValues[chan][samp][i]),
            print >> outFile        
        elif nuisTypes[i] == "gmN":
            print >> outFile, "{0:<14}{1:<6}{2:<4}".format(nuisLabels[i],nuisTypes[i],norms[i]),
            for chan in cardChannels:
                for samp in cardLabels:
                    if nuisValues[chan][samp][i] != "-" and float(nuisValues[chan][samp][i]) == 0 :
                        print >> outFile, "{0:^8.3f}".format(getYield(f,cardMassSamples[mass][samp],['ww'+chan],-1,"CONVLHT")/float(norms[i])),
                    else:
                        print >> outFile, "{0:^8}".format(nuisValues[chan][samp][i]),
            print >> outFile        
            
    
