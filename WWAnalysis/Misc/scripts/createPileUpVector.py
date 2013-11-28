#!/usr/bin/env python

from ROOT import TH1F,TFile
import os

# 2011 A
f = TFile( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Scales/pu2011A.root" )
h = f.Get("pileup")
dataCertified = [ h.GetBinContent(i)/h.GetSumOfWeights() for i in range(1,h.GetNbinsX()+1) ]
f.Close()

fout = open( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/python/certifiedPileUp_cfi.py",'w' )

print >> fout, "import FWCore.ParameterSet.Config as cms"
print >> fout
print >> fout, "pu2011A = [",
for data in dataCertified:
    print >> fout, data, ", ",
print >> fout, "]"
print >> fout

# 2011 B
f = TFile( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Scales/pu2011B.root" )
h = f.Get("pileup")
dataCertified = [ h.GetBinContent(i)/h.GetSumOfWeights() for i in range(1,h.GetNbinsX()+1) ]
f.Close()


print >> fout
print >> fout, "pu2011B = [",
for data in dataCertified:
    print >> fout, data, ", ",
print >> fout, "]"
print >> fout

# 2011 A+B
f = TFile( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Scales/pu2011AB.root" )
h = f.Get("pileup")
dataCertified = [ h.GetBinContent(i)/h.GetSumOfWeights() for i in range(1,h.GetNbinsX()+1) ]
f.Close()


print >> fout
print >> fout, "pu2011AB = [",
for data in dataCertified:
    print >> fout, data, ", ",
print >> fout, "]"
print >> fout

# MC
f = TFile( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Scales/s4MCPileUp.root" )
h = f.Get("hNPU")
dataCertified = [ h.GetBinContent(i)/h.GetSumOfWeights() for i in range(1,h.GetNbinsX()+1) ]
f.Close()


print >> fout
print >> fout, "puS4fromMC = [",
for data in dataCertified:
    print >> fout, data, ", ",
print >> fout, "]"
print >> fout

f = TFile( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Scales/s6MCPileUp.root" )
h = f.Get("pileup")
dataCertified = [ h.GetBinContent(i)/h.GetSumOfWeights() for i in range(1,h.GetNbinsX()+1) ]
f.Close()


print >> fout
print >> fout, "puS6fromMC = [",
for data in dataCertified:
    print >> fout, data, ", ",
print >> fout, "]"
print >> fout

fout.close()
