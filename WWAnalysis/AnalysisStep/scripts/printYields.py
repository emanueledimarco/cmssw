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

f = TFile(sys.argv[1])
yieldDir = f.Get('eventHists/yields')

# channels = ['wwelel']
# types = ['0']
channels = ['wwmumu','wwmuel','wwelmu','wwelel']
# types = ['0', 'IDCBL','ISOCBL','CONVCBL','IPCBL',
#               'IDCBT','ISOCBT','CONVCBT','IPCBT',
#               'IDLHL','ISOLHL','CONVLHL','IPLHL',
#               'IDLHT','ISOLHT','CONVLHT','IPLHT']
# types = ['IDLHT','ISOLHT','CONVLHT','IPLHT','ISOPFLHT','CONVPFLHT','IPPFLHT']
# types = ['0', 'IDLHT','ISOLHT','CONVLHT','IPLHT']
# types = ['0']
# types = ['CONVLHT']
types = ['1Jet']
nbins = 20

# hists = [ "004.DYtoMuMuS3", "101160.ggToH160toWWto2L2Nu", "023.TTJetsMad" ]
hists = [ "101160_ggToH160toWWto2L2Nu" ]

for hist in hists:
    for type in types:
        # print seperator
        print "{0:*^100}".format(type+" "+hist)
    
        # print header
        print "{0:^40}".format("Cut"),
        for chan in channels:
            print "{0:^16}".format(chan),
        print "\n{0:-^100}".format("")
    
        # print yields
        for bin in range(1,nbins+1):
            for chan in channels:
                localDir = yieldDir.Get(chan+type)
                h = localDir.Get(hist)
                if chan == channels[0]:
                    print "{0:^40}".format(h.GetXaxis().GetBinLabel(bin)),
                print "{0:^16.0f}".format(h.GetBinContent(bin)),
            print
    
        # print seperator
        print "{0:*^100}\n".format("")

# print "Yield Breakdowns @ {0} pb^-1".format(lumi)
# print "Channel,Data,Data Error,SM Signal,SM Error,4th Gen Signal,4th Gen Error,Background,Background Error"
# for channel in channels:
#     print "{0},{1},{2},{3},{4}".format(
#         channel,
#         myFormat(dataYields[channel][nbins],math.sqrt(dataErr2[channel][nbins])),
#         myFormat(sigYields[channel][nbins], math.sqrt(sigErr2[channel][nbins])),
#         myFormat(enYields[channel][nbins],  math.sqrt(enErr2[channel][nbins])),
#         myFormat(bkgYields[channel][nbins], math.sqrt(bkgErr2[channel][nbins]))
#     )
# print
# 
# for channel in channels:
#     print "Cut Breakdown by number of events for {0}.".format(channel)
#     for line in linesRaw[channel]:
#         print linesRaw[channel][line]
#     print
#     
#     print "Cut Breakdown by yield to {0} pb^-1 for {1}.".format(lumi,channel)
#     for line in linesScaled[channel]:
#         print linesScaled[channel][line]
#     print
#  
